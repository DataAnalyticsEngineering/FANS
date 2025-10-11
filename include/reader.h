
#ifndef READER_H
#define READER_H

#include <map>
#include <string>
#include <vector>
#include "mixedBCs.h"

using namespace std;

class Reader {
  public:
    // Default constructor
    Reader();

    // Destructor to free allocated memory
    ~Reader();

    // contents of input file:
    char             ms_filename[4096];    // Name of Micro-structure hdf5 file
    char             ms_datasetname[4096]; // Absolute path of Micro-structure in hdf5 file
    char             results_prefix[4096];
    int              n_mat;
    json             materialProperties;
    double           TOL;
    json             errorParameters;
    json             microstructure;
    int              n_it;
    vector<LoadCase> load_cases;
    string           problemType;
    string           matmodel;
    string           method;

    vector<string> resultsToWrite;

    // contents of microstructure file:
    vector<int>     dims;
    vector<double>  l_e;
    vector<double>  L;
    unsigned short *ms; // Micro-structure
    bool            is_zyx = true;

    int world_rank;
    int world_size;

    ptrdiff_t alloc_local;
    ptrdiff_t local_n0;
    ptrdiff_t local_0_start; // this is the x-value of the start point, not the index in the array
    ptrdiff_t local_n1;
    ptrdiff_t local_1_start;

    // void Setup(ptrdiff_t howmany);
    void ReadInputFile(char fn[]);
    void ReadMS(int hm);
    void ComputeVolumeFractions();
    // void ReadHDF5(char file_name[], char dset_name[]);
    void safe_create_group(hid_t file, const char *const name);

    template <typename T>
    void WriteSlab(T *data, int _howmany, const char *file_name, const char *dset_name);

    template <typename T>
    void WriteData(T *data, const char *file_name, const char *dset_name, hsize_t *dims, int rank);
};

template <typename T>
void Reader::WriteData(T *data, const char *file_name, const char *dset_name, hsize_t *dims, int rank)
{
    hid_t data_type;
    if (std::is_same<T, double>::value) {
        data_type = H5T_NATIVE_DOUBLE;
    } else if (std::is_same<T, unsigned char>::value) {
        data_type = H5T_NATIVE_UCHAR;
    } else {
        throw std::invalid_argument("Conversion of this data type to H5 data type not yet implemented");
    }

    hid_t plist_id;
    hid_t file_id;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    // TODO: refactor this into a general error handling method
    /* Save old error handler */
    herr_t (*old_func)(hid_t, void *);
    void *old_client_data;
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
    /* Turn off error handling */
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    file_id = H5Fopen(file_name, H5F_ACC_RDWR, plist_id);
    /* Restore previous error handler */
    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

    if (file_id < 0) {
        file_id = H5Fcreate(file_name, H5F_ACC_EXCL, H5P_DEFAULT, plist_id);
    }
    H5Pclose(plist_id);

    // Ensure all groups in the path are created
    safe_create_group(file_id, dset_name);

    // Check if the dataset already exists
    if (H5Lexists(file_id, dset_name, H5P_DEFAULT) > 0) {
        // If it exists, delete it
        H5Ldelete(file_id, dset_name, H5P_DEFAULT);
    }

    // Create the data space for the dataset
    hid_t dataspace_id = H5Screate_simple(rank, dims, NULL);
    if (dataspace_id < 0) {
        H5Fclose(file_id);
        throw std::runtime_error("Error creating dataspace");
    }

    // Create the dataset with default properties
    hid_t dataset_id = H5Dcreate2(file_id, dset_name, data_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id < 0) {
        H5Sclose(dataspace_id);
        H5Fclose(file_id);
        throw std::runtime_error("Error creating dataset");
    }

    // Write the data to the dataset
    if (H5Dwrite(dataset_id, data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0) {
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
        H5Fclose(file_id);
        throw std::runtime_error("Error writing data to dataset");
    }

    // Close the dataset and the file
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Fclose(file_id);
}

// this function has to be here because of the template
/* ---------------------------------------------------------------------------
 * Write a 4-D slab (local_n0 × Ny × Nz × howmany) from THIS MPI rank into an
 * HDF5 dataset whose global layout on disk is
 *
 *            Z  Y  X  k   (i.e. "zyx" + extra dim k)
 *
 * The caller supplies the data in logical order
 *            X  Y  Z  k .
 *
 * We therefore transpose in-memory once per write call.
 * --------------------------------------------------------------------------*/
template <typename T>
void Reader::WriteSlab(
    T          *data,     // in:  local slab, layout [X][Y][Z][k]
    int         _howmany, // global size of the 4th axis (k)
    const char *file_name,
    const char *dset_name)
{
    /*------------------------------------------------------------------*/
    /* 0. map C++ type -> native HDF5 type                              */
    /*------------------------------------------------------------------*/
    hid_t data_type;
    if (std::is_same<T, double>())
        data_type = H5T_NATIVE_DOUBLE;
    else if (std::is_same<T, float>())
        data_type = H5T_NATIVE_FLOAT;
    else if (std::is_same<T, unsigned char>())
        data_type = H5T_NATIVE_UCHAR;
    else if (std::is_same<T, int>())
        data_type = H5T_NATIVE_INT;
    else if (std::is_same<T, unsigned short>())
        data_type = H5T_NATIVE_USHORT;
    else
        throw std::invalid_argument("WriteSlab: unsupported data type");

    /*------------------------------------------------------------------*/
    /* 1. open or create the HDF5 file                                  */
    /*------------------------------------------------------------------*/
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    // H5Pset_fapl_mpio(plist_id, comm, info);   /* if you need MPI-IO */

    /* temporarily silence “file not found” during H5Fopen */
    herr_t (*old_func)(hid_t, void *);
    void *old_client_data;
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
    H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);

    hid_t file_id = H5Fopen(file_name, H5F_ACC_RDWR, plist_id);
    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data); /* restore */
    if (file_id < 0)                                     /* create if absent */
        file_id = H5Fcreate(file_name, H5F_ACC_EXCL, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    /*------------------------------------------------------------------*/
    /* 2. create the dataset (global dims =  Z Y X k ) if necessary     */
    /*------------------------------------------------------------------*/
    const hsize_t Nx   = static_cast<hsize_t>(dims[0]);
    const hsize_t Ny   = static_cast<hsize_t>(dims[1]);
    const hsize_t Nz   = static_cast<hsize_t>(dims[2]);
    const hsize_t kDim = static_cast<hsize_t>(_howmany);

    const int rank        = 4;
    hsize_t   dimsf[rank] = {Nz, Ny, Nx, kDim}; /* <<< Z Y X k */

    safe_create_group(file_id, dset_name);

    /* silence “dataset not found” during open */
    H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
    H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);

    hid_t dset_id = H5Dopen2(file_id, dset_name, H5P_DEFAULT);
    H5Eset_auto(H5E_DEFAULT, old_func, old_client_data); /* restore */

    if (dset_id < 0) {
        hid_t filespace = H5Screate_simple(rank, dimsf, nullptr);
        dset_id         = H5Dcreate2(file_id, dset_name, data_type,
                                     filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);

        /*--------------------------------------------------------------*/
        /*  add the attribute  permute_order = "zyx"                    */
        /*--------------------------------------------------------------*/
        const char perm_str[] = "zyx";
        hid_t      atype      = H5Tcopy(H5T_C_S1);
        H5Tset_size(atype, 4); /* 3 chars + '\0' */
        hid_t aspace = H5Screate(H5S_SCALAR);
        hid_t attr   = H5Acreate2(dset_id, "permute_order",
                                  atype, aspace, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, atype, perm_str);
        H5Aclose(attr);
        H5Sclose(aspace);
        H5Tclose(atype);
    }

    /*------------------------------------------------------------------*/
    /* 3. build FILE hyperslab  (slice along file-dim 2 = X axis)       */
    /*------------------------------------------------------------------*/
    hsize_t fcount[4]  = {Nz, Ny, static_cast<hsize_t>(local_n0), kDim};
    hsize_t foffset[4] = {0, 0, static_cast<hsize_t>(local_0_start), 0};

    hid_t filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, foffset, nullptr, fcount, nullptr);

    /*------------------------------------------------------------------*/
    /* 4. transpose local slab  [X][Y][Z][k]  ->  [Z][Y][X][k]          */
    /*------------------------------------------------------------------*/
    const size_t NxLoc = static_cast<size_t>(local_n0);
    const size_t NyLoc = static_cast<size_t>(Ny);
    const size_t NzLoc = static_cast<size_t>(Nz);
    const size_t kLoc  = static_cast<size_t>(kDim);

    size_t         slabElems = NxLoc * NyLoc * NzLoc * kLoc;
    std::vector<T> tmp(slabElems); /* automatic RAII buffer */

    for (size_t x = 0; x < NxLoc; ++x)
        for (size_t y = 0; y < NyLoc; ++y)
            for (size_t z = 0; z < NzLoc; ++z) {
                size_t srcBase = (((x * NyLoc) + y) * NzLoc + z) * kLoc; /* X-major */
                size_t dstBase = (((z * NyLoc) + y) * NxLoc + x) * kLoc; /* Z-major */
                std::memcpy(&tmp[dstBase], &data[srcBase], kLoc * sizeof(T));
            }

    /*------------------------------------------------------------------*/
    /* 5. MEMORY dataspace matches FILE slab exactly                    */
    /*------------------------------------------------------------------*/
    hid_t memspace = H5Screate_simple(4, fcount, nullptr);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    // H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    herr_t status = H5Dwrite(dset_id, data_type,
                             memspace, filespace, plist_id, tmp.data());
    if (status < 0)
        throw std::runtime_error("WriteSlab: H5Dwrite failed");

    /*------------------------------------------------------------------*/
    /* 6. tidy up                                                       */
    /*------------------------------------------------------------------*/
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Fclose(file_id);
}

#endif
