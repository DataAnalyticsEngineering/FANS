#include "solver.h"

template <int howmany, int n_str>
using MacroVector = Matrix<double, n_str, 1>;

template <int howmany, int n_str>
using ElementVector = Matrix<double, howmany * 8, 1>;

template <int howmany, int n_str, int padding>
static void tangent_residual(Solver<howmany, n_str> &solver, RealArray &residual, RealArray &delta_u,
                             const double *tangent, const MacroVector<howmany, n_str> &macro)
{
    auto *model = solver.matmanager->models[0];
    solver.template compute_residual_basic<padding>(
        residual, delta_u, [&](ElementVector<howmany, n_str> &ue, int, ptrdiff_t element) -> ElementVector<howmany, n_str> & {
            return model->element_tangent_residual(ue, macro, tangent + element * model->n_gp * n_str * n_str);
        });
}

template <int howmany, int n_str>
static MacroVector<howmany, n_str> tangent_average(Solver<howmany, n_str> &solver, RealArray &delta_u, RealArray &scratch,
                                                   const double *tangent, const MacroVector<howmany, n_str> &macro)
{
    auto                         *model   = solver.matmanager->models[0];
    MacroVector<howmany, n_str>   average = MacroVector<howmany, n_str>::Zero();
    ElementVector<howmany, n_str> zero    = ElementVector<howmany, n_str>::Zero();
    solver.template compute_residual_basic<0>(scratch, delta_u, [&](ElementVector<howmany, n_str> &ue, int, ptrdiff_t element) -> ElementVector<howmany, n_str> & {
        average += model->element_tangent_stress(ue, macro, tangent + element * model->n_gp * n_str * n_str);
        return zero;
    });
    MPI_Allreduce(MPI_IN_PLACE, average.data(), n_str, MPI_DOUBLE, MPI_SUM, solver.communicator);
    return average / (solver.n_x * solver.n_y * solver.n_z);
}

template <int howmany, int n_str>
static void solve_tangent(Solver<howmany, n_str> &solver, RealArray &delta_u, RealArray &r, RealArray &d,
                          RealArray &Kd, RealArray &z, const double *tangent, const MacroVector<howmany, n_str> &macro)
{
    tangent_residual<howmany, n_str, 2>(solver, r, delta_u, tangent, macro);
    d.setZero();
    double       rho = 0.0, error = std::sqrt(solver.dotProduct(r, r));
    const double limit = solver.reader.errorParameters["type"] == "relative" ? solver.TOL * error : solver.TOL;
    for (int i = 0; i < solver.n_it && error > limit; ++i) {
        solver.apply_preconditioner(r.data(), z.data());
        z *= -1.0;
        const double old_rho = rho;
        rho                  = solver.dotProduct(r, z);
        d                    = z + (i ? rho / old_rho : 0.0) * d;
        tangent_residual<howmany, n_str, 0>(solver, Kd, d, tangent, MacroVector<howmany, n_str>::Zero());
        const double alpha = rho / solver.dotProduct(d, Kd);
        delta_u -= alpha * d;
        r -= alpha * Kd;
        error = std::sqrt(solver.dotProduct(r, r));
    }
}

template <int howmany, int n_str>
MatrixXd compute_consistent_homogenized_tangent(Solver<howmany, n_str> &solver)
{
    if (solver.isMixedBCActive())
        throw std::runtime_error("Homogenized tangent is not supported with mixed boundary conditions.");
    const auto     n_elem = solver.local_n0 * solver.n_y * solver.n_z;
    const size_t   n_work = (solver.local_n0 + 1) * solver.n_y * solver.n_z * howmany;
    const size_t   n_fft  = std::max<size_t>(solver.reader.alloc_local * 2, (solver.local_n0 + 1) * solver.n_y * (solver.n_z + 2) * howmany);
    vector<double> displacement(n_work), search_direction(n_work);
    unique_ptr<double, decltype(&fftw_free)> residual(fftw_alloc_real(n_fft), fftw_free), product(fftw_alloc_real(n_fft), fftw_free);
    vector<double> local_tangent;
    auto          *model = solver.matmanager->models[0];
    model->compute_tangent_field(n_elem, solver.ms, local_tangent);
    RealArray delta_u(displacement.data(), solver.n_z * howmany, solver.local_n0 * solver.n_y, OuterStride<>(solver.n_z * howmany));
    RealArray r(residual.get(), solver.n_z * howmany, solver.local_n0 * solver.n_y, OuterStride<>((solver.n_z + 2) * howmany));
    RealArray d(search_direction.data(), solver.n_z * howmany, solver.local_n0 * solver.n_y, OuterStride<>(solver.n_z * howmany));
    RealArray Kd(product.get(), solver.n_z * howmany, solver.local_n0 * solver.n_y, OuterStride<>(solver.n_z * howmany));
    RealArray z(product.get(), solver.n_z * howmany, solver.local_n0 * solver.n_y, OuterStride<>((solver.n_z + 2) * howmany));
    MatrixXd  homogenized_tangent(n_str, n_str);
    for (int j = 0; j < n_str; ++j) {
        delta_u.setZero();
        const MacroVector<howmany, n_str> direction = MacroVector<howmany, n_str>::Unit(j);
        solve_tangent(solver, delta_u, r, d, Kd, z, local_tangent.data(), direction);
        homogenized_tangent.col(j) = tangent_average(solver, delta_u, Kd, local_tangent.data(), direction);
    }
    return homogenized_tangent;
}

template MatrixXd compute_consistent_homogenized_tangent(Solver<1, 3> &);
template MatrixXd compute_consistent_homogenized_tangent(Solver<3, 6> &);
template MatrixXd compute_consistent_homogenized_tangent(Solver<3, 9> &);
