#ifndef SOLVER_CG_H
#define SOLVER_CG_H

#include "solver.h"

template <int howmany>
class SolverCG : public Solver<howmany>{
public:
    using Solver<howmany>::n_x;
    using Solver<howmany>::n_y;
    using Solver<howmany>::n_z;
    using Solver<howmany>::local_n0;
    using Solver<howmany>::local_n1;
    using Solver<howmany>::v_u_real;
    using Solver<howmany>::v_r_real;

    SolverCG(Reader reader, Matmodel<howmany>* matmodel);

    double *s;
    double *d;
    double *rnew;
    RealArray s_real;
    RealArray d_real;
    RealArray rnew_real;

    void internalSolve();
    void LineSearchSecant();
    double dotProduct(RealArray& a, RealArray& b);

protected:
    using Solver<howmany>::iter;
};


template<int howmany>
SolverCG<howmany> :: SolverCG(Reader reader, Matmodel<howmany>* mat) : Solver<howmany>(reader, mat),

    s(fftw_alloc_real(reader.alloc_local * 2)),
    s_real(s, n_z * howmany, local_n0 * n_y, OuterStride<>((n_z + 2) * howmany)),

    rnew(fftw_alloc_real((local_n0 + 1) * n_y * n_z * howmany)),
    rnew_real(rnew, n_z * howmany, local_n0 * n_y, OuterStride<>(n_z * howmany)),

    d(fftw_alloc_real((local_n0 + 1) * n_y * n_z * howmany)),
    d_real(d, n_z * howmany, local_n0 * n_y, OuterStride<>(n_z * howmany))
{
    this->CreateFFTWPlans(this->v_r, (fftw_complex*) s, s);
}


template<int howmany>
double SolverCG<howmany> :: dotProduct(RealArray& a, RealArray& b){
    double local_value = (a * b).sum();
    double result;
    MPI_Allreduce(&local_value, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return result;
}


template<int howmany>
void SolverCG<howmany> :: internalSolve (){
    if(this->world_rank == 0)
        printf("\n# Start FANS - Conjugate Gradient Solver \n");

    LinearModel<howmany>* linearModel = dynamic_cast<LinearModel<howmany>*>(this->matmodel);
    bool islinear = (linearModel == NULL) ? false : true;

    s_real.setZero();
    v_u_real.setZero();
    d_real.setZero();
    for(ptrdiff_t i = local_n0 * n_y * n_z * howmany; i < (local_n0 + 1) * n_y * n_z * howmany; i++){
        this->v_u[i] = 0;
        d[i] = 0;
    }

    this->template compute_residual<2>(v_r_real, d_real);
    
    iter = 0;
    double err_rel = this->compute_error(v_r_real);

    double delta, delta0, deltamid;
    delta = 1.0L;

    while ((iter < this->n_it) && (err_rel > this->TOL)){
        
        deltamid = dotProduct(v_r_real, s_real);

        this->convolution();

        s_real *= -1;
        delta0 = delta;
        delta = dotProduct(v_r_real, s_real);

        d_real = s_real + fmax(0, (delta-deltamid)/delta0) * d_real;

        if (islinear) {
            Matrix<double, howmany*8, 1> res_e;
            this->template compute_residual_basic<0>(rnew_real, d_real, 
            [&](Matrix<double, howmany*8, 1>& ue, int mat_index) -> Matrix<double, howmany*8, 1>& {
                res_e.noalias() = linearModel->phase_stiffness[mat_index] * ue;
                return res_e;
            });
            
            double alpha = delta / dotProduct(d_real, rnew_real);
            v_r_real -= alpha * rnew_real;
            v_u_real -= alpha * d_real;
        }
        else { LineSearchSecant(); }

        iter ++;
        err_rel = this->compute_error(v_r_real);
    }
    if(this->world_rank == 0)
        printf("# Complete FANS - Conjugate Gradient Solver \n");
}

template<int howmany>
void SolverCG<howmany> :: LineSearchSecant() {
    double err = 10.0; int MaxIter = 10; double tol = 1e-2; int _iter = 0;
    double alpha_new = 0.0001;
    double alpha_old = 0;

    double r1pd;
    double rpd = dotProduct(v_r_real, d_real);

    while (((_iter < MaxIter) && (err > tol)) || (_iter < 2)){

        v_u_real += d_real * (alpha_new - alpha_old);
        this->template compute_residual<0>(rnew_real, v_u_real);
        r1pd = dotProduct(rnew_real, d_real);

        alpha_old = alpha_new;
        alpha_new *= rpd/(rpd-r1pd);

        err = fabs(alpha_new - alpha_old);
        _iter ++;
    }
    v_u_real += d_real * (alpha_new - alpha_old);
    v_r_real = rnew_real;
    if(this->world_rank == 0) printf("line search iter %i, alpha %f - error %e - ", _iter, alpha_new, err);
}
#endif