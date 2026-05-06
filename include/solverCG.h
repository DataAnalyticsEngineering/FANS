#ifndef SOLVER_CG_H
#define SOLVER_CG_H

#include "solver.h"
#include "logging.h"

template <int howmany, int n_str>
class SolverCG : public Solver<howmany, n_str> {
  public:
    using Solver<howmany, n_str>::n_x;
    using Solver<howmany, n_str>::n_y;
    using Solver<howmany, n_str>::n_z;
    using Solver<howmany, n_str>::local_n0;
    using Solver<howmany, n_str>::local_n1;
    using Solver<howmany, n_str>::v_u_real;
    using Solver<howmany, n_str>::v_r_real;

    SolverCG(Reader &reader, MaterialManager<howmany, n_str> *matmanager);
    virtual ~SolverCG();

    double   *s;
    double   *d;
    double   *rnew;
    RealArray s_real;
    RealArray d_real;
    RealArray rnew_real;
    double    alpha_warm{0.1};
    bool      ls_converged{true};

    void   internalSolve();
    void   LineSearchSecant();
    double dotProduct(RealArray &a, RealArray &b);

  protected:
    using Solver<howmany, n_str>::iter;
};

template <int howmany, int n_str>
SolverCG<howmany, n_str>::SolverCG(Reader &reader, MaterialManager<howmany, n_str> *matmgr)
    : Solver<howmany, n_str>(reader, matmgr),

      s(fftw_alloc_real(reader.alloc_local * 2)),
      s_real(s, n_z * howmany, local_n0 * n_y, OuterStride<>((n_z + 2) * howmany)),

      rnew(fftw_alloc_real((local_n0 + 1) * n_y * n_z * howmany)),
      rnew_real(rnew, n_z * howmany, local_n0 * n_y, OuterStride<>(n_z * howmany)),

      d(fftw_alloc_real((local_n0 + 1) * n_y * n_z * howmany)),
      d_real(d, n_z * howmany, local_n0 * n_y, OuterStride<>(n_z * howmany))
{
    this->CreateFFTWPlans(this->v_r, (fftw_complex *) s, s);
}

template <int howmany, int n_str>
double SolverCG<howmany, n_str>::dotProduct(RealArray &a, RealArray &b)
{
    double local_value = (a * b).sum();
    double result;
    MPI_Allreduce(&local_value, &result, 1, MPI_DOUBLE, MPI_SUM, this->communicator);
    return result;
}

template <int howmany, int n_str>
void SolverCG<howmany, n_str>::internalSolve()
{
    Log::solver->info() << "\t# Start FANS - Conjugate Gradient Solver \n";

    bool islinear = this->matmanager->all_linear;
    alpha_warm    = 0.1;

    s_real.setZero();
    d_real.setZero();
    for (ptrdiff_t i = local_n0 * n_y * n_z * howmany; i < (local_n0 + 1) * n_y * n_z * howmany; i++) {
        d[i] = 0;
    }

    this->template compute_residual<2>(v_r_real, v_u_real);

    iter           = 0;
    double err_rel = this->compute_error(v_r_real);

    double delta, delta0, deltamid;
    delta = 1.0L;

    while ((iter < this->n_it) && (err_rel > this->TOL)) {

        deltamid = dotProduct(v_r_real, s_real);

        this->convolution();

        s_real *= -1;
        delta0 = delta;
        delta  = dotProduct(v_r_real, s_real);

        if (islinear && !this->isMixedBCActive()) {
            d_real = s_real + fmax(0.0, (delta - deltamid) / delta0) * d_real;
            Matrix<double, howmany * 8, 1> res_e;
            this->template compute_residual_basic<0>(rnew_real, d_real,
                                                     [&](Matrix<double, howmany * 8, 1> &ue, int phase_id, ptrdiff_t element_idx) -> Matrix<double, howmany * 8, 1> & {
                                                         const MaterialInfo<howmany, n_str> &info = this->matmanager->get_info(phase_id);
                                                         res_e.noalias()                          = info.linear_model->phase_stiffness[info.local_mat_id] * ue;
                                                         return res_e;
                                                     });

            double alpha = delta / dotProduct(d_real, rnew_real);
            v_r_real -= alpha * rnew_real;
            v_u_real -= alpha * d_real;
        } else {
            d_real = s_real + (ls_converged ? fmax(0.0, (delta - deltamid) / delta0) : 0.0) * d_real;
            LineSearchSecant();
        }

        iter++;
        err_rel = this->compute_error(v_r_real);

        if (iter >= 2 && this->err_all[iter] > this->err_all[iter - 1] && this->err_all[iter - 1] > this->err_all[iter - 2])
            ls_converged = false; // Force a CG restart
    }
    Log::solver->info() << "# Complete FANS - Conjugate Gradient Solver \n";
}

template <int howmany, int n_str>
void SolverCG<howmany, n_str>::LineSearchSecant()
{
    const int    MaxIter    = this->reader.ls_max_iter;
    const double tol        = this->reader.ls_tol;
    int          _iter      = 0;
    double       alpha_prev = 0.0;
    double       alpha_curr = alpha_warm;

    double rpd = dotProduct(v_r_real, d_real);
    if (rpd >= 0.0) {
        ls_converged = false;
        alpha_warm   = 0.1;
        return;
    }
    v_u_real += d_real * alpha_curr;
    this->updateMixedBC();
    this->template compute_residual<0>(rnew_real, v_u_real);
    double       r1pd  = dotProduct(rnew_real, d_real);
    const double r1pd0 = this->isMixedBCActive() ? r1pd : rpd;

    double denom, alpha_next;
    while (_iter < MaxIter && fabs(r1pd) > tol * fabs(r1pd0)) {
        denom = r1pd - rpd;
        if (fabs(denom) < 1e-14 * (fabs(r1pd) + fabs(rpd)))
            break;

        alpha_next = alpha_curr - r1pd * (alpha_curr - alpha_prev) / denom;
        if (alpha_next <= 0.0 || alpha_next > 10.0)
            alpha_next = 0.5 * (alpha_prev + alpha_curr);

        v_u_real += d_real * (alpha_next - alpha_curr);
        alpha_prev = alpha_curr;
        rpd        = r1pd;
        alpha_curr = alpha_next;
        _iter++;

        this->template compute_residual<0>(rnew_real, v_u_real);
        r1pd = dotProduct(rnew_real, d_real);
    }
    ls_converged = (fabs(r1pd) <= tol * fabs(r1pd0));
    if (ls_converged) {
        alpha_warm = alpha_curr;
        v_r_real   = rnew_real;
    } else {
        v_u_real -= d_real * alpha_curr;
        alpha_warm = 0.1;
    }

    Log::solver->info() << Log::format("line search iter %i, alpha %f - error %e - ", _iter, alpha_curr, fabs(r1pd0) > 0.0 ? fabs(r1pd) / fabs(r1pd0) : 0.0);
}

template <int howmany, int n_str>
SolverCG<howmany, n_str>::~SolverCG()
{
    if (s) {
        fftw_free(s);
        s = nullptr;
    }
    if (rnew) {
        fftw_free(rnew);
        rnew = nullptr;
    }
    if (d) {
        fftw_free(d);
        d = nullptr;
    }
}
#endif
