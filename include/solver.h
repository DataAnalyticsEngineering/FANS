#ifndef SOLVER_H
#define SOLVER_H

#include "matmodel.h"

typedef Map<Array<double, Dynamic, Dynamic>, Unaligned, OuterStride<>> RealArray;

template <int howmany>
class Solver {
  public:
    Solver(Reader reader, Matmodel<howmany> *matmodel);

    Reader reader;

    const int world_rank;
    const int world_size;

    const ptrdiff_t n_x, n_y, n_z;
    // NOTE: the order in the declaration is very important because it is the same order in which the later initialization via member initializer lists takes place
    //  see https://stackoverflow.com/questions/1242830/constructor-initialization-list-evaluation-order
    const ptrdiff_t local_n0;
    const ptrdiff_t local_0_start; // this is the x-index of the start point, not the index in the array
    const ptrdiff_t local_n1;
    const ptrdiff_t local_1_start;

    const int          n_it;     //!< Max number of FANS iterations
    const double       TOL;      //!< Tolerance on relative error norm
    Matmodel<howmany> *matmodel; //!< Material Model

    unsigned char *ms;  // Micro-structure Binary
    double        *v_r; //!< Residual vector
    double        *v_u;
    double        *buffer_padding;

    RealArray      v_r_real; // can't do the "classname()" intialization here, and Map doesn't have a default constructor
    RealArray      v_u_real;
    Map<VectorXcd> rhat;

    ArrayXd                          err_all; //!< Absolute error history
    Matrix<double, howmany, Dynamic> fundamentalSolution;

    template<int padding, typename F>
    void iterateCubes(F f);

    void         solve();
    virtual void internalSolve(){}; // important to have "{}" here, otherwise we get an error about undefined reference to vtable

    template <int padding, typename F>
    void compute_residual_basic(RealArray &r_matrix, RealArray &u_matrix, F f);
    template <int padding>
    void compute_residual(RealArray &r_matrix, RealArray &u_matrix);

    void postprocess(Reader reader, const char resultsFileName[], int load_idx, int time_idx); //!< Computes Strain and stress

    void   convolution();
    double compute_error(RealArray &r);
    void   CreateFFTWPlans(double *in, fftw_complex *transformed, double *out);

    vector<double> homogenized_stress;

    // Homogenized stress which will be accessed from the library
    vector<double> get_homogenized_stress();

  protected:
    fftw_plan planfft, planifft;
    clock_t   fft_time, buftime;
    size_t    iter;
};

template <int howmany>
Solver<howmany>::Solver(Reader reader, Matmodel<howmany> *mat)
    : reader(reader),
      matmodel(mat),
      world_rank(reader.world_rank),
      world_size(reader.world_size),
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

      rhat((std::complex<double> *) v_r, local_n1 * n_x * (n_z / 2 + 1) * howmany), // actual initialization is below
      buffer_padding(fftw_alloc_real(n_y * (n_z + 2) * howmany))
{
    v_u_real.setZero();
    for (ptrdiff_t i = local_n0 * n_y * n_z * howmany; i < (local_n0 + 1) * n_y * n_z * howmany; i++) {
        this->v_u[i] = 0;
    }

    matmodel->initializeInternalVariables(local_n0 * n_y * n_z, 8);

    if (world_rank == 0) {
        printf("\n# Start creating Fundamental Solution(s) \n");
    }
    clock_t tot_time = clock();

    Matrix<double, howmany * 8, howmany * 8> Ker0 = matmodel->Compute_Reference_ElementStiffness();

    complex<double> tpi = 2 * acos(-1) * complex<double>(0, 1); //=2*pi*i

    auto etax = [&](double i) { return exp(tpi * i / (double) n_x); };
    auto etay = [&](double i) { return exp(tpi * i / (double) n_y); };
    auto etaz = [&](double i) { return exp(tpi * i / (double) n_z); };

    Matrix<complex<double>, 8, 1>    A;
    Matrix<double, 8, 8>             AA;
    Matrix<double, howmany, howmany> block;
    fundamentalSolution = Matrix<double, howmany, Dynamic>(howmany, (local_n1 * n_x * (n_z / 2 + 1) * (howmany + 1)) / 2);

    for (int i_y = 0; i_y < local_n1; ++i_y) {
        for (int i_x = 0; i_x < n_x; ++i_x) {
            for (int i_z = 0; i_z < n_z / 2 + 1; ++i_z) {
                if (i_x == 0 && local_1_start + i_y == 0 && i_z == 0) {
                    fundamentalSolution.template middleCols<howmany>(0).template triangularView<Lower>() = MatrixXd::Zero(howmany, howmany).triangularView<Lower>();
                    continue;
                }
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
                    fundamentalSolution.template middleCols<howmany>((ind / 2) * (howmany + 1)).template triangularView<Lower>() = block.inverse().template triangularView<Lower>();
                } else {
                    fundamentalSolution.template middleCols<howmany>((ind / 2) * (howmany + 1) + 1).template triangularView<Upper>() = block.inverse().template triangularView<Upper>();
                }
            }
        }
    }
    // Divided by n_el to scale the Fundamental solution so explicit normalization is not needed for FFT and IFFT
    fundamentalSolution /= (double) (n_x * n_y * n_z);

    tot_time = clock() - tot_time;
    if (world_rank == 0) {
        printf("# Complete; Time for construction of Fundamental Solution(s): %f seconds\n", double(tot_time) / CLOCKS_PER_SEC);
    }
}

template <int howmany>
void Solver<howmany>::CreateFFTWPlans(double *in, fftw_complex *transformed, double *out)
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
    planfft              = fftw_mpi_plan_many_dft_r2c(rank, n, howmany, iblock, oblock, in, transformed, MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);
    planifft             = fftw_mpi_plan_many_dft_c2r(rank, n, howmany, iblock, oblock, transformed, out, MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);

    // see https://eigen.tuxfamily.org/dox/group__TutorialMapClass.html#title3
    new (&rhat) Map<VectorXcd>((std::complex<double> *) transformed, local_n1 * n_x * (n_z / 2 + 1) * howmany);
}

// TODO: possibly circumvent the padding problem by accessing r as a matrix?
template <int howmany>
template <int padding, typename F>
void Solver<howmany>::compute_residual_basic(RealArray &r_matrix, RealArray &u_matrix, F f)
{

    double *r = r_matrix.data();
    double *u = u_matrix.data();
    r_matrix.setZero();
    // TODO: define another eigen Map for setting this part to zero?
    for (ptrdiff_t i = local_n0 * n_y * (n_z + padding) * howmany; i < (local_n0 + 1) * n_y * (n_z + padding) * howmany; i++) {
        r[i] = 0;
    }

    // int MPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf,
    //           int recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status)
    MPI_Sendrecv(u, n_y * n_z * howmany, MPI_DOUBLE, (world_rank + world_size - 1) % world_size, 0,
                 u + local_n0 * n_y * n_z * howmany, n_y * n_z * howmany, MPI_DOUBLE, (world_rank + 1) % world_size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

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
                 buffer_padding, n_y * (n_z + padding) * howmany, MPI_DOUBLE, (world_rank + world_size - 1) % world_size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    RealArray b(buffer_padding, n_z * howmany, n_y, OuterStride<>((n_z + padding) * howmany)); // NOTE: for any padding of more than 2, the buffer_padding has to be extended

    r_matrix.block(0, 0, n_z * howmany, n_y) += b; // matrix.block(i,j,p,q); is the block of size (p,q), starting at (i,j)
}

template <int howmany>
template <int padding>
void Solver<howmany>::compute_residual(RealArray &r_matrix, RealArray &u_matrix)
{
    compute_residual_basic<padding>(r_matrix, u_matrix, [&](Matrix<double, howmany * 8, 1> &ue, int mat_index, ptrdiff_t element_idx) -> Matrix<double, howmany * 8, 1> & {
        return matmodel->element_residual(ue, mat_index, element_idx);
    });
}

template <int howmany>
void Solver<howmany>::solve()
{

    err_all          = ArrayXd::Zero(n_it + 1);
    fft_time         = 0.0;
    clock_t tot_time = clock();
    internalSolve();
    tot_time = clock() - tot_time;
    // if( VERBOSITY > 5 ){
    if (world_rank == 0) {
        printf("# FFT Time per iteration .......   %2.6f sec\n", double(fft_time) / CLOCKS_PER_SEC / iter);
        printf("# Total FFT Time ...............   %2.6f sec\n", double(fft_time) / CLOCKS_PER_SEC);
        printf("# Total Time per iteration .....   %2.6f sec\n", double(tot_time) / CLOCKS_PER_SEC / iter);
        printf("# Total Time ...................   %2.6f sec\n", double(tot_time) / CLOCKS_PER_SEC);
        printf("# FFT contribution to total time   %2.6f %% \n", 100. * double(fft_time) / double(tot_time));
    }
    matmodel->updateInternalVariables();
}

template <int howmany>
template <int padding, typename F>
void Solver<howmany>::iterateCubes(F f)
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

template <int howmany>
void Solver<howmany>::convolution()
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

template <int howmany>
double Solver<howmany>::compute_error(RealArray &r)
{
    double err, err0, err_local;

    err_local = r.matrix().lpNorm<Infinity>();
    MPI_Allreduce(&err_local, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    err_all[iter]  = err;
    err0           = err_all[0];
    double err_rel = (iter == 0 ? 100 : err / err0);

    if (world_rank == 0) {
        if (iter == 0) {
            printf("Before 1st iteration: %16.8e\n", err0);
        } else if (iter == 1) {
            printf("it %3lu .... err %16.8e  / %8.4e, ratio: ------------- , FFT time: %2.6f sec\n", iter, err, err / err0, double(buftime) / CLOCKS_PER_SEC);
        } else {
            printf("it %3lu .... err %16.8e  / %8.4e, ratio: %4.8e, FFT time: %2.6f sec\n", iter, err, err / err0, err / err_all[iter - 1], double(buftime) / CLOCKS_PER_SEC);
        }
    }

    return err; // returns absolute error
    // return err_rel; // returns relative error
}

template <int howmany>
void Solver<howmany>::postprocess(Reader reader, const char resultsFileName[], int load_idx, int time_idx)
{
    int      n_str          = matmodel->n_str;
    VectorXd strain         = VectorXd::Zero(local_n0 * n_y * n_z * n_str);
    VectorXd stress         = VectorXd::Zero(local_n0 * n_y * n_z * n_str);
    VectorXd stress_average = VectorXd::Zero(n_str);
    VectorXd strain_average = VectorXd::Zero(n_str);

    // Initialize per-phase accumulators
    int              n_mat = reader.n_mat;
    vector<VectorXd> phase_stress_average(n_mat, VectorXd::Zero(n_str));
    vector<VectorXd> phase_strain_average(n_mat, VectorXd::Zero(n_str));
    vector<int>      phase_counts(n_mat, 0);

    MPI_Sendrecv(v_u, n_y * n_z * howmany, MPI_DOUBLE, (world_rank + world_size - 1) % world_size, 0,
                 v_u + local_n0 * n_y * n_z * howmany, n_y * n_z * howmany, MPI_DOUBLE, (world_rank + 1) % world_size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    Matrix<double, howmany * 8, 1> ue;
    int                            mat_index;
    iterateCubes<0>([&](ptrdiff_t *idx, ptrdiff_t *idxPadding) {
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < howmany; ++j) {
                ue(howmany * i + j, 0) = v_u[howmany * idx[i] + j];
            }
        }
        mat_index = ms[idx[0]];

        matmodel->getStrainStress(strain.segment(n_str * idx[0], n_str).data(), stress.segment(n_str * idx[0], n_str).data(), ue, mat_index, idx[0]);
        stress_average += stress.segment(n_str * idx[0], n_str);
        strain_average += strain.segment(n_str * idx[0], n_str);

        phase_stress_average[mat_index] += stress.segment(n_str * idx[0], n_str);
        phase_strain_average[mat_index] += strain.segment(n_str * idx[0], n_str);
        phase_counts[mat_index]++;
    });

    MPI_Allreduce(MPI_IN_PLACE, stress_average.data(), n_str, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, strain_average.data(), n_str, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    stress_average /= (n_x * n_y * n_z);
    strain_average /= (n_x * n_y * n_z);

    // Reduce per-phase accumulations across all processes
    for (int mat_index = 0; mat_index < n_mat; ++mat_index) {
        MPI_Allreduce(MPI_IN_PLACE, phase_stress_average[mat_index].data(), n_str, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, phase_strain_average[mat_index].data(), n_str, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &phase_counts[mat_index], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // Compute average for each phase
        if (phase_counts[mat_index] > 0) {
            phase_stress_average[mat_index] /= phase_counts[mat_index];
            phase_strain_average[mat_index] /= phase_counts[mat_index];
        }
    }

    if (world_rank == 0) {
        printf("# Effective Stress .. (");
        for (int i = 0; i < n_str; ++i)
            printf("%+f ", stress_average[i]);
        printf(") \n");
        printf("# Effective Strain .. (");
        for (int i = 0; i < n_str; ++i)
            printf("%+f ", strain_average[i]);
        printf(") \n\n");
    }

    // Write results to results h5 file
    auto writeData = [&](const char *resultName, const char *resultPrefix, auto *data, hsize_t size) {
        if (std::find(reader.resultsToWrite.begin(), reader.resultsToWrite.end(), resultName) != reader.resultsToWrite.end()) {
            char name[5096];
            sprintf(name, "%s/load%i/time_step%i/%s", reader.ms_datasetname, load_idx, time_idx, resultPrefix);
            hsize_t dims[1] = {size};
            reader.WriteData(data, resultsFileName, name, dims, 1);
        }
    };

    auto writeSlab = [&](const char *resultName, const char *resultPrefix, auto *data, int size) {
        if (std::find(reader.resultsToWrite.begin(), reader.resultsToWrite.end(), resultName) != reader.resultsToWrite.end()) {
            char name[5096];
            sprintf(name, "%s/load%i/time_step%i/%s", reader.ms_datasetname, load_idx, time_idx, resultPrefix);
            reader.WriteSlab(data, size, resultsFileName, name);
        }
    };

    for (int i = 0; i < world_size; ++i) {
        if (i == world_rank) {
            if (world_rank == 0) {
                writeData("stress_average", "stress_average", stress_average.data(), n_str);
                writeData("strain_average", "strain_average", strain_average.data(), n_str);
                writeData("absolute_error", "absolute_error", err_all.data(), iter + 1);

                for (int mat_index = 0; mat_index < n_mat; ++mat_index) {
                    char stress_name[512];
                    char strain_name[512];
                    sprintf(stress_name, "phase_stress_average_phase%d", mat_index);
                    sprintf(strain_name, "phase_strain_average_phase%d", mat_index);
                    writeData("phase_stress_average", stress_name, phase_stress_average[mat_index].data(), n_str);
                    writeData("phase_strain_average", strain_name, phase_strain_average[mat_index].data(), n_str);
                }
            }
            writeSlab("microstructure", "microstructure", ms, 1);
            writeSlab("displacement", "displacement", v_u, howmany);
            writeSlab("residual", "residual", v_r, howmany);
            writeSlab("strain", "strain", strain.data(), n_str);
            writeSlab("stress", "stress", stress.data(), n_str);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    matmodel->postprocess(*this, reader, resultsFileName, load_idx, time_idx);

    // Copy computed average strain to member variable
    homogenized_stress = stress_average;
}

template <int howmany>
std::vector<double> get_homogenized_stress()
{
    return homogenous_stress;
}

#endif
