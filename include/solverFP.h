#ifndef SOLVER_FP_H
#define SOLVER_FP_H

#include "solver.h"

template <int howmany, int n_str>
class SolverFP : public Solver<howmany, n_str> {
  public:
    using Solver<howmany, n_str>::n_x;
    using Solver<howmany, n_str>::n_y;
    using Solver<howmany, n_str>::n_z;
    using Solver<howmany, n_str>::local_n0;
    using Solver<howmany, n_str>::local_n1;
    using Solver<howmany, n_str>::v_u_real;
    using Solver<howmany, n_str>::v_r_real;

    SolverFP(Reader &reader, MaterialManager<howmany, n_str> *matmanager);

    void internalSolve();

  protected:
    using Solver<howmany, n_str>::iter;
};

template <int howmany, int n_str>
SolverFP<howmany, n_str>::SolverFP(Reader &reader, MaterialManager<howmany, n_str> *matmanager)
    : Solver<howmany, n_str>(reader, matmanager)
{
    this->CreateFFTWPlans(this->v_r, (fftw_complex *) this->v_r, this->v_r);
}

template <int howmany, int n_str>
void SolverFP<howmany, n_str>::internalSolve()
{
    if (this->world_rank == 0)
        printf("\n# Start FANS - Fixed Point Solver \n");

    this->template compute_residual<2>(v_r_real, v_u_real);

    iter           = 0;
    double err_rel = this->compute_error(v_r_real);

    while ((iter < this->n_it) && (err_rel > this->TOL)) {

        this->convolution();
        v_u_real -= v_r_real;
        this->updateMixedBC();
        this->template compute_residual<2>(v_r_real, v_u_real);

        iter++;
        err_rel = this->compute_error(v_r_real);
    }
    if (this->world_rank == 0)
        printf("# Complete FANS - Fixed Point Solver \n");
}
#endif
