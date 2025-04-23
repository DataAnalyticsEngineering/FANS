#include "general.h"
#include "reader.h"

#include "H5Cpp.h"
#include "fftw3-mpi.h"
#include "hdf5.h"
#include "mpi.h"
#include "stdlib.h"

#include "H5FDmpi.h"
#include "H5FDmpio.h"

void Reader::ComputeVolumeFractions()
{
    int    local_max  = INT_MIN;
    int    local_min  = INT_MAX;
    size_t local_size = local_n0 * dims[1] * dims[2];

    // Find the local maximum and minimum material indices
    for (size_t i = 0; i < local_size; i++) {
        if (ms[i] > local_max) {
            local_max = ms[i];
        }
        if (ms[i] < local_min) {
            local_min = ms[i];
        }
    }

    // Find the global maximum and minimum material indices
    int global_max, global_min;
    MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    // Calculate total number of materials (accounting for negative indices)
    n_mat = global_max - global_min + 1;

    if (world_rank == 0) {
        printf("# Number of materials: %i (from %i to %i)\n", n_mat, global_min, global_max);
        printf("# Volume fractions\n");
    }

    // Using dynamic allocation for arrays since we don't know size at compile time
    std::vector<long>   vol_frac(n_mat, 0);
    std::vector<double> v_frac(n_mat, 0.0);

    for (size_t i = 0; i < local_size; i++) {
        int index = ms[i] - global_min; // Adjust index to start from 0
        vol_frac[index]++;
    }

    for (int i = 0; i < n_mat; i++) {
        long vf;
        MPI_Allreduce(&(vol_frac[i]), &vf, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
        v_frac[i] = double(vf) / double(dims[0] * dims[1] * dims[2]);
        if (world_rank == 0)
            printf("# material %4i    vol. frac. %10.4f%%  \n", i + global_min, 100. * v_frac[i]);
    }
}

void Reader ::ReadInputFile(char fn[])
{
    try {

        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        ifstream i(fn);
        json     j;
        i >> j;

        microstructure = j["microstructure"];
        strcpy(ms_filename, microstructure["filepath"].get<string>().c_str());
        strcpy(ms_datasetname, microstructure["datasetname"].get<string>().c_str());
        L = microstructure["L"].get<vector<double>>();

        if (j.contains("results_prefix")) {
            strcpy(results_prefix, j["results_prefix"].get<string>().c_str());
        } else {
            strcpy(results_prefix, "");
        }

        errorParameters = j["error_parameters"];
        TOL             = errorParameters["tolerance"].get<double>();
        n_it            = j["n_it"].get<int>();
        g0              = j["macroscale_loading"].get<vector<vector<vector<double>>>>();

        problemType = j["problem_type"].get<string>();
        matmodel    = j["matmodel"].get<string>();
        method      = j["method"].get<string>();

        json j_mat     = j["material_properties"];
        resultsToWrite = j["results"].get<vector<string>>(); // Read the results_to_write field

        if (world_rank == 0) {
            printf("# microstructure file name: \t '%s'\n", ms_filename);
            printf("# microstructure dataset name: \t '%s'\n", ms_datasetname);
            printf(
                "# FANS error measure: \t %s %s error  \n",
                errorParameters["type"].get<string>().c_str(),
                errorParameters["measure"].get<string>().c_str());
            printf("# FANS Tolerance: \t %10.5e\n", errorParameters["tolerance"].get<double>());
            printf("# Max iterations: \t %6i\n", n_it);
        }

        for (auto it = j_mat.begin(); it != j_mat.end(); ++it) {
            materialProperties[it.key()] = it.value();

            if (world_rank == 0) {
                cout << "# " << it.key() << ":\t ";
                if (it.value().is_array()) {
                    for (const auto &elem : it.value()) {
                        if (elem.is_number()) {
                            printf("   %10.5f", elem.get<double>());
                        } else if (elem.is_string()) {
                            cout << "   " << elem.get<string>();
                        } else {
                            cout << "   " << elem;
                        }
                    }
                } else if (it.value().is_number()) {
                    printf("   %10.5f", it.value().get<double>());
                } else if (it.value().is_string()) {
                    cout << "   " << it.value().get<string>();
                } else {
                    cout << "   " << it.value();
                }
                printf("\n");
            }
        }

    } catch (const std::exception &e) {
        fprintf(stderr, "ERROR trying to read input file '%s' for FANS\n", fn);
        exit(10);
    }
}

void Reader::safe_create_group(hid_t file, const char *const name)
{
    // no leading '/' --> exit
    const char DELIMITER = '/';
    if (name[0] != DELIMITER)
        return;

    // copy name to buffer
    char buffer[4096];
    strcpy(buffer, name);
    char *str = buffer;
    str       = strchr(str + 1, DELIMITER);
    while (str != NULL) {
        // while another / character is found
        long int l = str - buffer; // length of substring
        buffer[l]  = '\0';         // temporary 'end of string'

        // safely create the group if needed
        hid_t group;

        /* Save old error handler */
        herr_t (*old_func)(hid_t, void *);
        void *old_client_data;
        H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
        /* Turn off error handling */
        H5Eset_auto(H5E_DEFAULT, NULL, NULL);

        group = H5Gopen(file, buffer, H5P_DEFAULT);
        /* Restore previous error handler */
        H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

        if (group < 0) {
            group = H5Gcreate(file, buffer, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        }
        H5Gclose(group);

        buffer[l] = DELIMITER; // restore original string

        str = strchr(str + 1, DELIMITER); // find next delimiter
    }
}

void Reader ::ReadMS(int hm)
{

    hid_t   file_id, dset_id;    /* file and dataset identifiers */
    hid_t   filespace, memspace; /* file and memory dataspace identifiers */
    hid_t   data_type;
    hsize_t _dims[3]; /* dataset dimensions */
    hsize_t count[3]; /* hyperslab selection parameters */
    hsize_t offset[3];
    hid_t   plist_id; /* property list identifier */
    herr_t  status;

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Info info = MPI_INFO_NULL;

    // Set up file access property list with parallel I/O access
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    // H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);        // "set File Access Property List"

    // Open the file collectively and release property list identifier.
    file_id = H5Fopen(ms_filename, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);

    // Create property list for collective dataset write.
    // plist_id = H5Pcreate(H5P_DATASET_XFER);
    // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);   // "set Data Transfer Property List" (x means transfer)
    plist_id = H5P_DEFAULT;

    dset_id = H5Dopen2(file_id, ms_datasetname, plist_id);

    hid_t dspace = H5Dget_space(dset_id);
    int   rank   = H5Sget_simple_extent_dims(dspace, _dims, NULL);
    data_type    = H5T_NATIVE_INT; // H5Dget_type(dset_id);

    dims.resize(3);
    dims[0] = _dims[0];
    dims[1] = _dims[1];
    dims[2] = _dims[2];

    l_e.resize(3);
    l_e[0] = L[0] / double(dims[0]);
    l_e[1] = L[1] / double(dims[1]);
    l_e[2] = L[2] / double(dims[2]);

    if (world_rank == 0) {
        printf("# grid size set to [%i x %i x %i] --> %i voxels \nMicrostructure length: [%3.6f x %3.6f x %3.6f]\n", dims[0], dims[1], dims[2], dims[0] * dims[1] * dims[2], L[0], L[1], L[2]);
        // if(dims[0] % 2 != 0)	fprintf(stderr, "[ FANS3D_Grid ] WARNING: n_x is not a multiple of 2\n");
        // if(dims[1] % 2 != 0)	fprintf(stderr, "[ FANS3D_Grid ] WARNING: n_y is not a multiple of 2\n");
        // if(dims[2] % 2 != 0)	fprintf(stderr, "[ FANS3D_Grid ] WARNING: n_z is not a multiple of 2\n");
        if (dims[0] / 4 < world_size)
            throw std::runtime_error("[ FANS3D_Grid ] ERROR: Please decrease the number of processes or increase the grid size to ensure that each process has at least 4 boxels in the x direction.");
        printf("Voxel length: [%1.8f, %1.8f, %1.8f]\n", l_e[0], l_e[1], l_e[2]);
    }

    const ptrdiff_t n[3]   = {dims[0], dims[1], dims[2] / 2 + 1};
    ptrdiff_t       block0 = FFTW_MPI_DEFAULT_BLOCK;
    ptrdiff_t       block1 = FFTW_MPI_DEFAULT_BLOCK;

    // see https://fftw.org/doc/Basic-and-advanced-distribution-interfaces.html
    // and https://www.fftw.org/fftw3_doc/Transposed-distributions.html
    // on there it is recommended to use one of fftw's allocation functions "to ensure optimal alignment"

    /* there is no documentation for this method, so here is the signature from "fftw3-mpi.h"
    FFTW_EXTERN ptrdiff_t XM(local_size_many_transposed)	\
     (int rnk, const ptrdiff_t *n, ptrdiff_t howmany,		\
      ptrdiff_t block0, ptrdiff_t block1, MPI_Comm comm,	\
      ptrdiff_t *local_n0, ptrdiff_t *local_0_start,		\
      ptrdiff_t *local_n1, ptrdiff_t *local_1_start);		\
    */

    alloc_local = fftw_mpi_local_size_many_transposed(rank, n, hm, block0, block1, MPI_COMM_WORLD, &local_n0, &local_0_start, &local_n1, &local_1_start);

    if (local_n0 < 4)
        throw std::runtime_error("[ FANS3D_Grid ] ERROR: Number of voxels in x-direction is less than 4 in process " + to_string(world_rank));
    MPI_Barrier(MPI_COMM_WORLD);

    // Each process defines a dataset in memory which reads a hyperslab from the file
    count[0] = local_n0;
    count[1] = dims[1];
    count[2] = dims[2];
    memspace = H5Screate_simple(rank, count, NULL);

    // Select hyperslab in the file.
    offset[0] = local_0_start;
    offset[1] = 0;
    offset[2] = 0;
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    ms     = FANS_malloc<int>(count[0] * count[1] * count[2]);
    status = H5Dread(dset_id, data_type, memspace, filespace, plist_id, this->ms);

    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Fclose(file_id);
    // H5Tclose(data_type);

    this->ComputeVolumeFractions();
}

// The code above is based on this example: Hyperslab_by_row.c

// /*
//  *  This example writes data to the HDF5 file by rows.
//  *  Number of processes is assumed to be 1 or multiples of 2 (up to 8)
//  */

// #include "hdf5.h"
// #include "stdlib.h"

// #define H5FILE_NAME     "SDS_row.h5"
// #define DATASETNAME 	"IntArray"
// #define NX     8                      /* dataset dimensions */
// #define NY     5
// #define RANK   2

// int
// main (int argc, char **argv)
// {
//     /*
//      * HDF5 APIs definitions
//      */
//     hid_t       file_id, dset_id;         /* file and dataset identifiers */
//     hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
//     hsize_t     dimsf[2];                 /* dataset dimensions */
//     int         *data;                    /* pointer to data buffer to write */
//     hsize_t	count[2];	          /* hyperslab selection parameters */
//     hsize_t	offset[2];
//     hid_t	plist_id;                 /* property list identifier */
//     int         i;
//     herr_t	status;

//     /*
//      * MPI variables
//      */
//     int mpi_size, mpi_rank;
//     MPI_Comm comm  = MPI_COMM_WORLD;
//     MPI_Info info  = MPI_INFO_NULL;

//     /*
//      * Initialize MPI
//      */
//     MPI_Init(&argc, &argv);
//     MPI_Comm_size(comm, &mpi_size);
//     MPI_Comm_rank(comm, &mpi_rank);

//     /*
//      * Set up file access property list with parallel I/O access
//      */
//      plist_id = H5Pcreate(H5P_FILE_ACCESS);
//      H5Pset_fapl_mpio(plist_id, comm, info);

//     /*
//      * Create a new file collectively and release property list identifier.
//      */
//     file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
//     H5Pclose(plist_id);

//     /*
//      * Create the dataspace for the dataset.
//      */
//     dimsf[0] = NX;
//     dimsf[1] = NY;
//     filespace = H5Screate_simple(RANK, dimsf, NULL);

//     /*
//      * Create the dataset with default properties and close filespace.
//      */
//     dset_id = H5Dcreate(file_id, DATASETNAME, H5T_NATIVE_INT, filespace,
// 			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     H5Sclose(filespace);

//     /*
//      * Each process defines dataset in memory and writes it to the hyperslab
//      * in the file.
//      */
//     count[0] = dimsf[0]/mpi_size;
//     count[1] = dimsf[1];
//     offset[0] = mpi_rank * count[0];
//     offset[1] = 0;
//     memspace = H5Screate_simple(RANK, count, NULL);

//     /*
//      * Select hyperslab in the file.
//      */
//     filespace = H5Dget_space(dset_id);
//     H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

//     /*
//      * Initialize data buffer
//      */
//     data = (int *) malloc(sizeof(int)*count[0]*count[1]);
//     for (i=0; i < count[0]*count[1]; i++) {
//         data[i] = mpi_rank + 10;
//     }

//     /*
//      * Create property list for collective dataset write.
//      */
//     plist_id = H5Pcreate(H5P_DATASET_XFER);
//     H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

//     status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
// 		      plist_id, data);
//     free(data);

//     /*
//      * Close/release resources.
//      */
//     H5Dclose(dset_id);
//     H5Sclose(filespace);
//     H5Sclose(memspace);
//     H5Pclose(plist_id);
//     H5Fclose(file_id);

//     MPI_Finalize();

//     return 0;
// }
