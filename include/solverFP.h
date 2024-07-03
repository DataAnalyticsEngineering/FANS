#ifndef SOLVER_FP_H
#define SOLVER_FP_H

#include "solver.h"

template <int howmany>
class SolverFP : public Solver<howmany>{
public:
    using Solver<howmany>::n_x;
    using Solver<howmany>::n_y;
    using Solver<howmany>::n_z;
    using Solver<howmany>::local_n0;
    using Solver<howmany>::local_n1;
    using Solver<howmany>::v_u_real;
    using Solver<howmany>::v_r_real;

    SolverFP(Reader reader, Matmodel<howmany>* matmodel);

    void internalSolve();
    
protected:
    using Solver<howmany>::iter;
};

template<int howmany>
SolverFP<howmany> :: SolverFP(Reader reader, Matmodel<howmany>* matmodel) : Solver<howmany>(reader, matmodel)
{
    this->CreateFFTWPlans(this->v_r, (fftw_complex*) this->v_r, this->v_r);
}   

template<int howmany>
void SolverFP<howmany> :: internalSolve() {
    if (this->world_rank == 0)
        printf("\n# Start FANS - Fixed Point Solver \n");

    v_u_real.setZero();
    for(ptrdiff_t i = local_n0 * n_y * n_z * howmany; i < (local_n0 + 1) * n_y * n_z * howmany; i++){
        this->v_u[i] = 0;
    }
    this->template compute_residual<2>(v_r_real, v_u_real);

    iter = 0;
    double err_rel = this->compute_error(v_r_real);

    while ((iter < this->n_it) && (err_rel > this->TOL)){

        this->convolution();
        v_u_real -= v_r_real;
        this->template compute_residual<2>(v_r_real, v_u_real);

        iter ++;
        err_rel = this->compute_error(v_r_real);
    }
    if (this->world_rank == 0)
        printf("# Complete FANS - Fixed Point Solver \n");
}
#endif