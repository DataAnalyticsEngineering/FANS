#ifndef SOLVER_H
#define SOLVER_H

#include "matmodel.h"
#include "MaterialManager.h"

typedef Map<Array<double, Dynamic, Dynamic>, Unaligned, OuterStride<>> RealArray;

template <int howmany, int n_str>
class Solver : private MixedBCController<howmany> {
  public:
    Solver(Reader &reader, MaterialManager<howmany, n_str> *matmanager);
    virtual ~Solver();

    Reader &reader;

    const int world_rank;
    const int world_size;
    MPI_Comm  communicator;

    const ptrdiff_t n_x, n_y, n_z;
    // NOTE: the order in the declaration is very important because it is the same order in which the later initialization via member initializer lists takes place
    //  see https://stackoverflow.com/questions/1242830/constructor-initialization-list-evaluation-order
    const ptrdiff_t local_n0;
    const ptrdiff_t local_0_start; // this is the x-index of the start point, not the index in the array
    const ptrdiff_t local_n1;
    const ptrdiff_t local_1_start;

    const int                        n_it;       //!< Max number of FANS iterations
    double                           TOL;        //!< Tolerance on relative error norm
    MaterialManager<howmany, n_str> *matmanager; //!< Material Manager

    unsigned short *ms;  // Micro-structure
    double         *v_r; //!< Residual vector
    double         *v_u;
    double         *v_u_prev; //!< Previous displacement for extrapolation
    double         *buffer_padding;

    RealArray      v_r_real; // can't do the "classname()" intialization here, and Map doesn't have a default constructor
    RealArray      v_u_real;
    Map<VectorXcd> rhat;

    ArrayXd                          err_all; //!< Absolute error history
    Matrix<double, howmany, Dynamic> fundamentalSolution;

    vector<vector<ptrdiff_t>> batch_elems;     //!< elements of each batching model
    vector<double>            batch_ue;        //!< nodal displacements of every element
    vector<double>            batch_gp_stress; //!< Gauss-point stresses of every element

    template <int padding, typename F>
    void iterateCubes(F f);
    void update_ghost_layer(double *u); //!< the next rank's first layer of nodal values, behind this rank's own

    void         solve();
    void         extrapolateDisplacement(); //!< Linear extrapolation for next time step
    virtual void internalSolve() {};        // important to have "{}" here, otherwise we get an error about undefined reference to vtable

    template <int padding, typename F>
    void compute_residual_basic(RealArray &r_matrix, RealArray &u_matrix, F f);
    template <int padding>
    void compute_residual(RealArray &r_matrix, RealArray &u_matrix);

    void postprocess(Reader &reader, int load_idx, int time_idx); //!< Computes Strain and stress

    void   convolution();
    double compute_error(RealArray &r, const std::string &details = {});
    void   CreateFFTWPlans(double *in, fftw_complex *transformed, double *out);

    VectorXd homogenized_strain;
    VectorXd homogenized_stress;
    VectorXd get_homogenized_stress();
    void     evaluate_batched_stress(double *u); //!< in src/plugin_material.cpp

    MatrixXd homogenized_tangent;
    MatrixXd get_homogenized_tangent(double pert_param); //!< in src/tangent.cpp

    void enableMixedBC(const MixedBC &mbc, size_t step)
    {
        this->activate(*this, mbc, step);
    }
    void disableMixedBC()
    {
        this->mixed_active = false;
    }
    bool isMixedBCActive()
    {
        return this->mixed_active;
    }
    void updateMixedBC()
    {
        this->update(*this);
    }

  private:
    template <typename _Matrix_Type_>
    inline _Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double tolerance) const
    {
        Eigen::JacobiSVD<_Matrix_Type_> svd(a, Eigen::ComputeFullU | Eigen::ComputeFullV);
        return svd.matrixV() *
               (svd.singularValues().array() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() *
               svd.matrixU().adjoint();
    }
    void computeFundamentalSolution();

  protected:
    fftw_plan planfft, planifft;
    clock_t   fft_time, buftime;
    size_t    iter;
};

template <int howmany, int n_str>
Solver<howmany, n_str>::Solver(Reader &reader, MaterialManager<howmany, n_str> *matmgr)
    : reader(reader),
      matmanager(matmgr),
      world_rank(reader.world_rank),
      world_size(reader.world_size),
      communicator(reader.communicator),
      n_x(reader.dims[0]),
      n_y(reader.dims[1]),
      n_z(reader.dims[2]),
      local_n0(reader.local_n0),
      local_n1(reader.local_n1),
      local_0_start(reader.local_0_start),
      local_1_start(reader.local_1_start),

      n_it(reader.n_it),
      TOL(reader.TOL),
      ms(reader.ms),

      v_r(fftw_alloc_real(std::max(reader.alloc_local * 2, (local_n0 + 1) * n_y * (n_z + 2) * howmany))),
      v_r_real(v_r, n_z * howmany, local_n0 * n_y, OuterStride<>((n_z + 2) * howmany)),

      v_u(fftw_alloc_real((local_n0 + 1) * n_y * n_z * howmany)),
      v_u_real(v_u, n_z * howmany, local_n0 * n_y, OuterStride<>(n_z * howmany)),
      v_u_prev(fftw_alloc_real(local_n0 * n_y * n_z * howmany)),

      rhat((std::complex<double> *) v_r, local_n1 * n_x * (n_z / 2 + 1) * howmany), // actual initialization is below
      buffer_padding(fftw_alloc_real(n_y * (n_z + 2) * howmany))
{
    v_u_real.setZero();
    for (ptrdiff_t i = local_n0 * n_y * n_z * howmany; i < (local_n0 + 1) * n_y * n_z * howmany; i++) {
        this->v_u[i] = 0;
    }
    std::memset(v_u_prev, 0, local_n0 * n_y * n_z * howmany * sizeof(double));

    matmanager->initialize_internal_variables(local_n0 * n_y * n_z, matmanager->models[0]->n_gp);

    computeFundamentalSolution();
}

template <int howmany, int n_str>
void Solver<howmany, n_str>::computeFundamentalSolution()
{
    Log::logger().info("# Start creating Fundamental Solution(s) ");
    clock_t tot_time = clock();

    Matrix<double, howmany * 8, howmany * 8> Ker0 = matmanager->models[0]->Compute_Reference_ElementStiffness(matmanager->kapparef_mat);

    complex<double> tpi = 2 * acos(-1) * complex<double>(0, 1); //=2*pi*i

    auto etax = [&](double i) { return exp(tpi * i / (double) n_x); };
    auto etay = [&](double i) { return exp(tpi * i / (double) n_y); };
    auto etaz = [&](double i) { return exp(tpi * i / (double) n_z); };

    Matrix<complex<double>, 8, 1>    A;
    Matrix<double, 8, 8>             AA;
    Matrix<double, howmany, howmany> block;
    fundamentalSolution = Matrix<double, howmany, Dynamic>(howmany, (local_n1 * n_x * (n_z / 2 + 1) * (howmany + 1)) / 2);
    fundamentalSolution.setZero();

    for (int i_y = 0; i_y < local_n1; ++i_y) {
        for (int i_x = 0; i_x < n_x; ++i_x) {
            for (int i_z = 0; i_z < n_z / 2 + 1; ++i_z) {
                if (i_x != 0 || (local_1_start + i_y) != 0 || i_z != 0) {

                    A(0, 0) = 1.0;
                    A(1, 0) = etax(i_x);
                    A(2, 0) = etay(local_1_start + i_y);
                    A(3, 0) = etax(i_x) * etay(local_1_start + i_y);
                    A(4, 0) = etaz(i_z);
                    A(5, 0) = etax(i_x) * etaz(i_z);
                    A(6, 0) = etaz(i_z) * etay(local_1_start + i_y);
                    A(7, 0) = etax(i_x) * etay(local_1_start + i_y) * etaz(i_z);
                    AA      = A.real() * A.real().transpose() + A.imag() * A.imag().transpose();

                    for (int i = 0; i < howmany; i++) {
                        for (int j = i; j < howmany; j++) {
                            block(i, j) = (Ker0.template block<8, 8>(8 * i, 8 * j).array() * AA.array()).sum();
                            block(j, i) = block(i, j); // we'd like to avoid this, but block.selfadjointView<Upper>().inverse() does not work
                        }
                    }
                    ptrdiff_t ind = i_y * n_x * (n_z / 2 + 1) + i_x * (n_z / 2 + 1) + i_z;
                    if (ind % 2 == 0) {
                        fundamentalSolution.template middleCols<howmany>((ind / 2) * (howmany + 1)).template triangularView<Lower>() = pseudoInverse(block, 1e-14).template triangularView<Lower>();
                    } else {
                        fundamentalSolution.template middleCols<howmany>((ind / 2) * (howmany + 1) + 1).template triangularView<Upper>() = pseudoInverse(block, 1e-14).template triangularView<Upper>();
                    }
                }
            }
        }
    }
    // // Divided by n_el to scale the Fundamental solution so explicit normalization is not needed for FFT and IFFT
    fundamentalSolution /= (double) (n_x * n_y * n_z);

    tot_time = clock() - tot_time;
    Log::logger().info("# Complete; Time for construction of Fundamental Solution(s): {:.6f} seconds", double(tot_time) / CLOCKS_PER_SEC);
}

template <int howmany, int n_str>
void Solver<howmany, n_str>::CreateFFTWPlans(double *in, fftw_complex *transformed, double *out)
{
    int       rank   = 3;
    ptrdiff_t iblock = FFTW_MPI_DEFAULT_BLOCK;
    ptrdiff_t oblock = FFTW_MPI_DEFAULT_BLOCK;

    // see https://fftw.org/doc/MPI-Plan-Creation.html
    // NOTE: according to https://fftw.org/doc/Multi_002ddimensional-MPI-DFTs-of-Real-Data.html:
    // "As for the serial transforms, the sizes you pass to the ‘plan_dft_r2c’ and ‘plan_dft_c2r’ are the n0 × n1 × n2 × … × nd-1 dimensions of the real data"
    // "That is, you call the appropriate ‘local size’ function for the n0 × n1 × n2 × … × (nd-1/2 + 1) complex data"
    // so we need to use a different n than for fftw_mpi_local_size_many !
    // But, according to https://fftw.org/doc/MPI-Plan-Creation.html the BLOCK sizes must be the same:
    // "These must be the same block sizes as were passed to the corresponding ‘local_size’ function"
    const ptrdiff_t n[3] = {n_x, n_y, n_z};
    planfft              = fftw_mpi_plan_many_dft_r2c(rank, n, howmany, iblock, oblock, in, transformed, communicator, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);
    planifft             = fftw_mpi_plan_many_dft_c2r(rank, n, howmany, iblock, oblock, transformed, out, communicator, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);

    // see https://eigen.tuxfamily.org/dox/group__TutorialMapClass.html#title3
    new (&rhat) Map<VectorXcd>((std::complex<double> *) transformed, local_n1 * n_x * (n_z / 2 + 1) * howmany);
}

// The elements of a rank's last layer reach into the first layer of nodes of the next rank
template <int howmany, int n_str>
void Solver<howmany, n_str>::update_ghost_layer(double *u)
{
    const int layer = n_y * n_z * howmany;
    MPI_Sendrecv(u, layer, MPI_DOUBLE, (world_rank + world_size - 1) % world_size, 0,
                 u + local_n0 * layer, layer, MPI_DOUBLE, (world_rank + 1) % world_size, 0, communicator, MPI_STATUS_IGNORE);
}

// TODO: possibly circumvent the padding problem by accessing r as a matrix?
template <int howmany, int n_str>
template <int padding, typename F>
void Solver<howmany, n_str>::compute_residual_basic(RealArray &r_matrix, RealArray &u_matrix, F f)
{

    double *r = r_matrix.data();
    double *u = u_matrix.data();
    r_matrix.setZero();
    // TODO: define another eigen Map for setting this part to zero?
    for (ptrdiff_t i = local_n0 * n_y * (n_z + padding) * howmany; i < (local_n0 + 1) * n_y * (n_z + padding) * howmany; i++) {
        r[i] = 0;
    }

    update_ghost_layer(u);

    Matrix<double, howmany * 8, 1> ue;

    iterateCubes<padding>([&](ptrdiff_t *idx, ptrdiff_t *idxPadding) {
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < howmany; j++) {
                ue(howmany * i + j, 0) = u[howmany * idx[i] + j] - u[howmany * idx[0] + j];
            }
        }
        Matrix<double, howmany * 8, 1> &res_e = f(ue, ms[idx[0]], idx[0]);

        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < howmany; j++) {
                r[howmany * idxPadding[i] + j] += res_e(howmany * i + j, 0);
            }
        }
    });

    MPI_Sendrecv(r + local_n0 * n_y * (n_z + padding) * howmany, n_y * (n_z + padding) * howmany, MPI_DOUBLE, (world_rank + 1) % world_size, 0,
                 buffer_padding, n_y * (n_z + padding) * howmany, MPI_DOUBLE, (world_rank + world_size - 1) % world_size, 0, communicator, MPI_STATUS_IGNORE);

    RealArray b(buffer_padding, n_z * howmany, n_y, OuterStride<>((n_z + padding) * howmany)); // NOTE: for any padding of more than 2, the buffer_padding has to be extended

    r_matrix.block(0, 0, n_z * howmany, n_y) += b; // matrix.block(i,j,p,q); is the block of size (p,q), starting at (i,j)
}

template <int howmany, int n_str>
template <int padding>
void Solver<howmany, n_str>::compute_residual(RealArray &r_matrix, RealArray &u_matrix)
{
    if (matmanager->any_batched)
        evaluate_batched_stress(u_matrix.data());
    compute_residual_basic<padding>(r_matrix, u_matrix, [&](Matrix<double, howmany * 8, 1> &ue, int phase_id, ptrdiff_t element_idx) -> Matrix<double, howmany * 8, 1> & {
        const MaterialInfo<howmany, n_str> &info = matmanager->get_info(phase_id);
        return info.model->element_residual(ue, info.local_mat_id, element_idx);
    });
}

template <int howmany, int n_str>
void Solver<howmany, n_str>::solve()
{

    err_all          = ArrayXd::Zero(n_it + 1);
    fft_time         = 0.0;
    clock_t tot_time = clock();
    internalSolve();
    tot_time = clock() - tot_time;
    Log::logger().info("# FFT Time per iteration .......   {:.6f} sec", iter == 0 ? 0.0 : double(fft_time) / CLOCKS_PER_SEC / iter);
    Log::logger().info("# Total FFT Time ...............   {:.6f} sec", double(fft_time) / CLOCKS_PER_SEC);
    Log::logger().info("# Total Time per iteration .....   {:.6f} sec", iter == 0 ? 0.0 : double(tot_time) / CLOCKS_PER_SEC / iter);
    Log::logger().info("# Total Time ...................   {:.6f} sec", double(tot_time) / CLOCKS_PER_SEC);
    Log::logger().info("# FFT contribution to total time   {:.6f} %", 100. * double(fft_time) / double(tot_time));
    if (matmanager->any_batched)
        evaluate_batched_stress(v_u); // the converged step: its stresses for the output, its state to keep
    matmanager->update_internal_variables();
}

template <int howmany, int n_str>
void Solver<howmany, n_str>::extrapolateDisplacement()
{
    const size_t n = local_n0 * n_y * n_z * howmany;
    for (size_t i = 0; i < n; ++i) {
        const double delta = v_u[i] - v_u_prev[i];
        v_u_prev[i]        = v_u[i];
        v_u[i] += delta; // v_u = v_u + (v_u - v_u_prev)
    }
}

template <int howmany, int n_str>
template <int padding, typename F>
void Solver<howmany, n_str>::iterateCubes(F f)
{

    auto Idx = [&](int i_x, int i_y) {
        if (i_y >= n_y)
            i_y -= n_y;
        return (n_z) * (n_y * i_x + i_y);
    };
    auto IdxPadding = [&](int i_x, int i_y) {
        if (i_y >= n_y)
            i_y -= n_y;
        return (n_z + padding) * (n_y * i_x + i_y);
    };
    ptrdiff_t idx[8], idxPadding[8];

    for (int i_x = 0; i_x < local_n0; ++i_x) {
        for (int i_y = 0; i_y < n_y; ++i_y) {

            idx[0] = Idx(i_x, i_y);
            idx[1] = Idx(i_x + 1, i_y);
            idx[2] = Idx(i_x, i_y + 1);
            idx[3] = Idx(i_x + 1, i_y + 1);
            idx[4] = idx[0] + 1;
            idx[5] = idx[1] + 1;
            idx[6] = idx[2] + 1;
            idx[7] = idx[3] + 1;

            idxPadding[0] = IdxPadding(i_x, i_y);
            idxPadding[1] = IdxPadding(i_x + 1, i_y);
            idxPadding[2] = IdxPadding(i_x, i_y + 1);
            idxPadding[3] = IdxPadding(i_x + 1, i_y + 1);
            idxPadding[4] = idxPadding[0] + 1;
            idxPadding[5] = idxPadding[1] + 1;
            idxPadding[6] = idxPadding[2] + 1;
            idxPadding[7] = idxPadding[3] + 1;

            for (int i_z = 0; i_z < n_z - 1; ++i_z) {
                f(idx, idxPadding);
                idx[0]++;
                idx[1]++;
                idx[2]++;
                idx[3]++;
                idx[4]++;
                idx[5]++;
                idx[6]++;
                idx[7]++;

                idxPadding[0]++;
                idxPadding[1]++;
                idxPadding[2]++;
                idxPadding[3]++;
                idxPadding[4]++;
                idxPadding[5]++;
                idxPadding[6]++;
                idxPadding[7]++;
            }

            idx[4] -= n_z;
            idx[5] -= n_z;
            idx[6] -= n_z;
            idx[7] -= n_z;

            idxPadding[4] -= n_z;
            idxPadding[5] -= n_z;
            idxPadding[6] -= n_z;
            idxPadding[7] -= n_z;

            f(idx, idxPadding);
        }
    }
}

template <int howmany, int n_str>
void Solver<howmany, n_str>::convolution()
{

    // it is important that at least one of the dimensions n_x and n_z is divisible by two (or local_n1, but that can't be guaranteed from the outside)
    // discussion of real times complex: https://forum.kde.org/viewtopic.php?f=74&t=85678

    clock_t dtime = clock();
    fftw_execute(planfft);
    fft_time += clock() - dtime;
    buftime = clock() - dtime;

    Matrix<complex<double>, howmany, howmany> tmp;
    for (ptrdiff_t i = 0; i < (local_n1 * n_x * (n_z / 2 + 1)) / 2; i++) {

        tmp                                          = fundamentalSolution.template middleCols<howmany>(i * (howmany + 1)).template cast<complex<double>>();
        rhat.segment<howmany>(2 * i * howmany)       = tmp.template selfadjointView<Lower>() * rhat.segment<howmany>(2 * i * howmany);
        tmp                                          = fundamentalSolution.template middleCols<howmany>(i * (howmany + 1) + 1).template cast<complex<double>>();
        rhat.segment<howmany>((2 * i + 1) * howmany) = tmp.template selfadjointView<Upper>() * rhat.segment<howmany>((2 * i + 1) * howmany);
    }

    dtime = clock();
    fftw_execute(planifft);
    fft_time += clock() - dtime;
    buftime += clock() - dtime;
}

template <int howmany, int n_str>
double Solver<howmany, n_str>::compute_error(RealArray &r, const std::string &details)
{
    double             err_local;
    const std::string &measure = reader.errorParameters["measure"].get<std::string>();
    if (measure == "L1") {
        err_local = r.matrix().lpNorm<1>();
    } else if (measure == "L2") {
        err_local = r.matrix().lpNorm<2>();
    } else if (measure == "Linfinity") {
        err_local = r.matrix().lpNorm<Infinity>();
    } else {
        throw std::runtime_error("Unknown measure type: " + measure);
    }

    double err;
    MPI_Allreduce(&err_local, &err, 1, MPI_DOUBLE, MPI_MAX, communicator);

    err_all[iter]  = err;
    double err0    = err_all[0];
    double err_rel = (err0 == 0.0 ? 0.0 : err / err0);

    if (iter == 0)
        Log::logger().info("Before 1st iteration: {:16.8e}", err0);
    else
        Log::logger().info("it {:3} .... err {:16.8e}  / {:8.4e}, ratio: {:4.8e}, FFT time: {:.6f} sec{}", iter, err, err_rel, iter == 1 || err_all[iter - 1] == 0.0 ? 0.0 : err / err_all[iter - 1], double(buftime) / CLOCKS_PER_SEC, details);

    const std::string &error_type = reader.errorParameters["type"].get<std::string>();
    if (error_type == "absolute") {
        return err;
    } else if (error_type == "relative") {
        return err_rel;
    } else {
        throw std::runtime_error("Unknown error type: " + error_type);
    }
}

template <int howmany, int n_str>
void Solver<howmany, n_str>::postprocess(Reader &reader, int load_idx, int time_idx)
{
    int n_gp = matmanager->models[0]->n_gp;

    // Check what user requested
    auto &results        = reader.resultsToWrite;
    bool  need_stress    = std::find(results.begin(), results.end(), "stress") != results.end();
    bool  need_stress_gp = std::find(results.begin(), results.end(), "stress_gp") != results.end();
    bool  need_strain    = std::find(results.begin(), results.end(), "strain") != results.end();
    bool  need_strain_gp = std::find(results.begin(), results.end(), "strain_gp") != results.end();

    bool need_global_avg = std::find(results.begin(), results.end(), "stress_average") != results.end() ||
                           std::find(results.begin(), results.end(), "strain_average") != results.end();

    bool need_phase_avg = false;
    for (int mat_idx = 0; mat_idx < reader.n_mat; ++mat_idx) {
        char name[512];
        sprintf(name, "phase_stress_average_phase%d", mat_idx);
        if (std::find(results.begin(), results.end(), name) != results.end()) {
            need_phase_avg = true;
            break;
        }
    }

    // Determine if we need to compute stress/strain at all
    bool need_compute = need_stress || need_stress_gp || need_strain || need_strain_gp || need_global_avg || need_phase_avg;

    // The requested fields; the others stay empty
    VectorXd strain_elem(need_strain ? local_n0 * n_y * n_z * n_str : 0);
    VectorXd stress_elem(need_stress ? local_n0 * n_y * n_z * n_str : 0);
    VectorXd strain_gp(need_strain_gp ? local_n0 * n_y * n_z * n_gp * n_str : 0);
    VectorXd stress_gp(need_stress_gp ? local_n0 * n_y * n_z * n_gp * n_str : 0);

    VectorXd stress_average = VectorXd::Zero(n_str);
    VectorXd strain_average = VectorXd::Zero(n_str);

    // Initialize per-phase accumulators
    int              n_mat = reader.n_mat;
    vector<VectorXd> phase_stress_average(n_mat, VectorXd::Zero(n_str));
    vector<VectorXd> phase_strain_average(n_mat, VectorXd::Zero(n_str));
    vector<int>      phase_counts(n_mat, 0);

    update_ghost_layer(v_u);

    Matrix<double, howmany * 8, 1> ue;
    int                            phase_id;

    if (need_compute) {
        iterateCubes<0>([&](ptrdiff_t *idx, ptrdiff_t *idxPadding) {
            for (int i = 0; i < 8; ++i) {
                for (int j = 0; j < howmany; ++j) {
                    ue(howmany * i + j, 0) = v_u[howmany * idx[i] + j];
                }
            }
            phase_id = ms[idx[0]];

            const MaterialInfo<howmany, n_str> &info = matmanager->get_info(phase_id);

            // Temporary storage for element averages
            double elem_strain_avg[n_str];
            double elem_stress_avg[n_str];

            // Compute once - populates internal eps/sigma
            info.model->getStrainStress(elem_strain_avg, elem_stress_avg, ue, info.local_mat_id, idx[0]);

            // Store element averages if requested
            if (need_stress) {
                for (int c = 0; c < n_str; ++c) {
                    stress_elem(idx[0] * n_str + c) = elem_stress_avg[c];
                }
            }
            if (need_strain) {
                for (int c = 0; c < n_str; ++c) {
                    strain_elem(idx[0] * n_str + c) = elem_strain_avg[c];
                }
            }

            // Copy all GP data if requested
            if (need_stress_gp) {
                memcpy(&stress_gp(idx[0] * n_gp * n_str),
                       info.model->get_sigma_data(),
                       n_gp * n_str * sizeof(double));
            }
            if (need_strain_gp) {
                memcpy(&strain_gp(idx[0] * n_gp * n_str),
                       info.model->get_eps_data(),
                       n_gp * n_str * sizeof(double));
            }

            // Accumulate for global and phase averages
            if (need_global_avg || need_phase_avg) {
                for (int c = 0; c < n_str; ++c) {
                    stress_average(c) += elem_stress_avg[c];
                    strain_average(c) += elem_strain_avg[c];
                    phase_stress_average[phase_id](c) += elem_stress_avg[c];
                    phase_strain_average[phase_id](c) += elem_strain_avg[c];
                }
                phase_counts[phase_id]++;
            }
        });
    }

    MPI_Allreduce(MPI_IN_PLACE, stress_average.data(), n_str, MPI_DOUBLE, MPI_SUM, communicator);
    MPI_Allreduce(MPI_IN_PLACE, strain_average.data(), n_str, MPI_DOUBLE, MPI_SUM, communicator);
    stress_average /= (n_x * n_y * n_z);
    strain_average /= (n_x * n_y * n_z);

    // Reduce per-phase accumulations across all processes
    for (int mat_index = 0; mat_index < n_mat; ++mat_index) {
        MPI_Allreduce(MPI_IN_PLACE, phase_stress_average[mat_index].data(), n_str, MPI_DOUBLE, MPI_SUM, communicator);
        MPI_Allreduce(MPI_IN_PLACE, phase_strain_average[mat_index].data(), n_str, MPI_DOUBLE, MPI_SUM, communicator);
        MPI_Allreduce(MPI_IN_PLACE, &phase_counts[mat_index], 1, MPI_INT, MPI_SUM, communicator);

        // Compute average for each phase
        if (phase_counts[mat_index] > 0) {
            phase_stress_average[mat_index] /= phase_counts[mat_index];
            phase_strain_average[mat_index] /= phase_counts[mat_index];
        }
    }

    if (Log::logger().should_log(spdlog::level::info)) {
        std::ostringstream output;
        output << std::showpos << std::scientific << std::setprecision(12)
               << "# Effective Stress .. (" << stress_average.transpose()
               << " ) \n# Effective Strain .. (" << strain_average.transpose() << " ) \n";
        Log::logger().info("{}", output.str());
    }
    homogenized_stress = stress_average;
    homogenized_strain = strain_average;

    // u_total = u + G X, with G the macroscale gradient of the nodal values: the temperature
    // gradient, the strain (ml holds its Mandel components) or F - I
    const vector<double>      &ml  = matmanager->models[0]->macroscale_loading;
    constexpr double           rs2 = 0.7071067811865475; // 1.0 / std::sqrt(2.0)
    Matrix<double, howmany, 3> G;
    if constexpr (n_str == 3)
        G << ml[0], ml[1], ml[2];
    else if constexpr (n_str == 6)
        G << ml[0], ml[3] * rs2, ml[4] * rs2,
            ml[3] * rs2, ml[1], ml[5] * rs2,
            ml[4] * rs2, ml[5] * rs2, ml[2];
    else
        G << ml[0] - 1.0, ml[1], ml[2],
            ml[3], ml[4] - 1.0, ml[5],
            ml[6], ml[7], ml[8] - 1.0;

    VectorXd  u_total(local_n0 * n_y * n_z * howmany);
    ptrdiff_t n = 0;
    for (ptrdiff_t ix = 0; ix < local_n0; ++ix) {
        const double x = (local_0_start + ix) * reader.l_e[0] - reader.L[0] / 2.0;
        for (ptrdiff_t iy = 0; iy < n_y; ++iy) {
            const double y = iy * reader.l_e[1] - reader.L[1] / 2.0;
            for (ptrdiff_t iz = 0; iz < n_z; ++iz, ++n) {
                const double z = iz * reader.l_e[2] - reader.L[2] / 2.0;
                for (int i = 0; i < howmany; ++i)
                    u_total[howmany * n + i] = v_u[howmany * n + i] + (G(i, 0) * x + G(i, 1) * y + G(i, 2) * z);
            }
        }
    }

    hsize_t dims[1] = {static_cast<hsize_t>(n_str)};
    reader.writeData("stress_average", load_idx, time_idx, stress_average.data(), dims, 1);
    reader.writeData("strain_average", load_idx, time_idx, strain_average.data(), dims, 1);
    for (int mat_index = 0; mat_index < n_mat; ++mat_index) {
        char stress_name[512];
        char strain_name[512];
        sprintf(stress_name, "phase_stress_average_phase%d", mat_index);
        sprintf(strain_name, "phase_strain_average_phase%d", mat_index);
        reader.writeData(stress_name, load_idx, time_idx, phase_stress_average[mat_index].data(), dims, 1);
        reader.writeData(strain_name, load_idx, time_idx, phase_strain_average[mat_index].data(), dims, 1);
    }
    dims[0] = iter + 1;
    reader.writeData("absolute_error", load_idx, time_idx, err_all.data(), dims, 1);

    vector<int> rank_field(local_n0 * n_y * n_z, world_rank);
    reader.writeSlab("mpi_rank", load_idx, time_idx, rank_field.data(), {1});
    reader.writeSlab("microstructure", load_idx, time_idx, ms, {1});
    reader.writeSlab("displacement_fluctuation", load_idx, time_idx, v_u, {howmany});
    reader.writeSlab("displacement", load_idx, time_idx, u_total.data(), {howmany});
    reader.writeSlab("residual", load_idx, time_idx, v_r, {howmany});

    if (need_strain)
        reader.writeSlab("strain", load_idx, time_idx, strain_elem.data(), {n_str});
    if (need_stress)
        reader.writeSlab("stress", load_idx, time_idx, stress_elem.data(), {n_str});
    if (need_strain_gp)
        reader.writeSlab("strain_gp", load_idx, time_idx, strain_gp.data(), {n_gp, n_str});
    if (need_stress_gp)
        reader.writeSlab("stress_gp", load_idx, time_idx, stress_gp.data(), {n_gp, n_str});

    matmanager->postprocess(*this, reader, load_idx, time_idx);

    // Compute homogenized tangent only if requested
    if (find(reader.resultsToWrite.begin(), reader.resultsToWrite.end(), "homogenized_tangent") != reader.resultsToWrite.end()) {
        homogenized_tangent = get_homogenized_tangent(1e-6);
        hsize_t dims[2]     = {static_cast<hsize_t>(n_str), static_cast<hsize_t>(n_str)};
        if (Log::logger().should_log(spdlog::level::info)) {
            std::ostringstream output;
            output << "# Homogenized tangent: \n"
                   << std::setprecision(12) << homogenized_tangent << '\n';
            Log::logger().info("{}", output.str());
        }
        reader.writeData("homogenized_tangent", load_idx, time_idx, homogenized_tangent.data(), dims, 2);
    }
}

template <int howmany, int n_str>
VectorXd Solver<howmany, n_str>::get_homogenized_stress()
{
    homogenized_stress = VectorXd::Zero(n_str);

    update_ghost_layer(v_u);
    if (matmanager->any_batched)
        evaluate_batched_stress(v_u);

    Matrix<double, howmany * 8, 1> ue;
    Matrix<double, n_str, 1>       strain, stress; // of one element
    int                            phase_id;
    iterateCubes<0>([&](ptrdiff_t *idx, ptrdiff_t *idxPadding) {
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < howmany; ++j) {
                ue(howmany * i + j, 0) = v_u[howmany * idx[i] + j];
            }
        }
        phase_id = ms[idx[0]];

        const MaterialInfo<howmany, n_str> &info = matmanager->get_info(phase_id);
        info.model->getStrainStress(strain.data(), stress.data(), ue, info.local_mat_id, idx[0]);
        homogenized_stress += stress;
    });

    MPI_Allreduce(MPI_IN_PLACE, homogenized_stress.data(), n_str, MPI_DOUBLE, MPI_SUM, communicator);
    homogenized_stress /= (n_x * n_y * n_z);

    return homogenized_stress;
}

template <int howmany, int n_str>
Solver<howmany, n_str>::~Solver()
{
    if (v_r) {
        fftw_free(v_r);
        v_r = nullptr;
    }
    if (v_u) {
        fftw_free(v_u);
        v_u = nullptr;
    }
    if (v_u_prev) {
        fftw_free(v_u_prev);
        v_u_prev = nullptr;
    }
    if (buffer_padding) {
        fftw_free(buffer_padding);
        buffer_padding = nullptr;
    }
    if (planfft) {
        fftw_destroy_plan(planfft);
    }
    if (planifft) {
        fftw_destroy_plan(planifft);
    }
}

#endif
