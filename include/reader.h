
#ifndef READER_H
#define READER_H

#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <vector>
#include <unsupported/Eigen/CXX11/Tensor>
#include "mixedBCs.h"
#include "hdf5.h"
#include "H5FDmpio.h"
#include "mpi.h"

using namespace std;

class Reader {
  public:
    // Default constructor
    Reader() = default;
    Reader(const MPI_Comm &comm);

    // Destructor to free allocated memory
    ~Reader();

    // contents of input file:
    char             ms_filename[4096]{};    // Name of Micro-structure hdf5 file
    char             ms_datasetname[4096]{}; // Absolute path of Micro-structure in hdf5 file
    char             results_prefix[4096]{};
    int              n_mat;
    json             inputJson; // Complete input JSON (for MaterialManager)
    json             materialProperties;
    double           TOL;
    json             errorParameters;
    int              ls_max_iter{5};
    double           ls_tol{1e-2};
    bool             extrapolate_displacement{true};
    json             microstructure;
    int              n_it;
    vector<LoadCase> load_cases;
    string           problemType;
    string           matmodel;
    string           method;
    string           strain_type{"small"}; // "small" (default) or "large"
    string           FE_type;              // "HEX8" (default), "HEX8R", or "BBAR"
    vector<string>   resultsToWrite;
    char             results_filename[4096]{}; // Output HDF5 filename
    char             dataset_name[8192]{};     // Base path for results in HDF5 file
    hid_t            results_file_id = -1;     // Open HDF5 file handle for results

    // contents of microstructure file:
    vector<int>     dims;
    vector<double>  l_e;
    vector<double>  L;
    unsigned short *ms{nullptr}; // Micro-structure
    bool            is_zyx = true;

    int      world_rank;
    int      world_size;
    MPI_Comm communicator;

    ptrdiff_t alloc_local;
    ptrdiff_t local_n0;
    ptrdiff_t local_0_start; // this is the x-value of the start point, not the index in the array
    ptrdiff_t local_n1;
    ptrdiff_t local_1_start;

    // void Setup(ptrdiff_t howmany);
    void ReadInputFile(const std::string &input_fn);
    void ReadMS(int hm);
    void FreeMS();
    void ComputeVolumeFractions();
    // void ReadHDF5(char file_name[], char dset_name[]);
    void safe_create_group(hid_t file, const char *const name);
    void OpenResultsFile(const char *output_fn); // Open results file once
    void CloseResultsFile();                     // Explicitly close results file

    // Convenience methods to check if a result should be written and write it
    template <typename T>
    void writeData(const char *fieldName, int load_idx, int time_idx, const T *data, const hsize_t *shape, int rank);
    template <typename T>
    void writeSlab(const char *fieldName, int load_idx, int time_idx, const T *data, const std::vector<int> &extra_dims);

    template <typename T>
    void WriteData(const T *data, const char *dset_name, const hsize_t *shape, int rank);
    template <typename T>
    void WriteSlab(const T *data, const std::vector<int> &extra_dims, const char *dset_name);

    // Read a dataset of any HDF5 file, on every rank: all of it, its shape, or
    // this rank's slab of a per-voxel dataset on the grid
    template <typename T>
    static std::vector<T>       ReadData(const string &file, const string &dset_name);
    static std::vector<hsize_t> DataShape(const string &file, const string &dset_name);
    template <typename T>
    void ReadSlab(T *data, const std::vector<int> &extra_dims, const string &file, const string &dset_name) const;

    // The microstructure's group, e.g. "/sphere/32x32x32/" for "/sphere/32x32x32/ms",
    // where the datasets that go with it are
    string MSGroup() const;
};

// Opens a dataset for reading; the caller closes it and file_id
hid_t h5_open_dataset(const string &file, const string &dset_name, hid_t &file_id);

// Whether a dataset is laid out [Z][Y][X]... (its permute_order attribute is
// "zyx", or missing) rather than [X][Y][Z]...
bool h5_is_zyx(hid_t dset_id);

// The native HDF5 type of T
template <typename T>
hid_t h5_native_type()
{
    if constexpr (std::is_same_v<T, double>)
        return H5T_NATIVE_DOUBLE;
    else if constexpr (std::is_same_v<T, float>)
        return H5T_NATIVE_FLOAT;
    else if constexpr (std::is_same_v<T, int>)
        return H5T_NATIVE_INT;
    else if constexpr (std::is_same_v<T, unsigned short>)
        return H5T_NATIVE_USHORT;
    else
        static_assert(!sizeof(T), "no native HDF5 type for T");
}

// [A][B][C][extra] -> [C][B][A][extra]: [X][Y][Z]... <-> [Z][Y][X]..., its own inverse
template <typename T>
void swap_xz(const T *in, T *out, Eigen::Index a, Eigen::Index b, Eigen::Index c, Eigen::Index extra)
{
    Eigen::TensorMap<Eigen::Tensor<const T, 4, Eigen::RowMajor>> from(in, a, b, c, extra);
    Eigen::TensorMap<Eigen::Tensor<T, 4, Eigen::RowMajor>>       to(out, c, b, a, extra);
    to = from.shuffle(Eigen::array<Eigen::Index, 4>{2, 1, 0, 3});
}

// All ranks call; rank 0's data is written, replacing an existing dataset
template <typename T>
void Reader::WriteData(const T *data, const char *dset_name, const hsize_t *shape, int rank)
{
    if (results_file_id < 0)
        throw std::runtime_error("WriteData: results file is not open");
    safe_create_group(results_file_id, dset_name);
    if (H5Lexists(results_file_id, dset_name, H5P_DEFAULT) > 0)
        H5Ldelete(results_file_id, dset_name, H5P_DEFAULT);

    hid_t space   = H5Screate_simple(rank, shape, nullptr);
    hid_t dset_id = H5Dcreate2(results_file_id, dset_name, h5_native_type<T>(), space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (world_rank != 0)
        H5Sselect_none(space); // the others take part, writing nothing
    hid_t xfer = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(xfer, H5FD_MPIO_INDEPENDENT);
    const herr_t status = H5Dwrite(dset_id, h5_native_type<T>(), H5S_ALL, space, xfer, world_rank == 0 ? data : nullptr);
    H5Pclose(xfer);
    H5Dclose(dset_id);
    H5Sclose(space);
    if (status < 0)
        throw std::runtime_error(string("WriteData: cannot write ") + dset_name);
}

// All ranks call with their slab, data as [X][Y][Z][d1][d2]... (local_n0 x Ny x Nz x ...);
// on disk [Z][Y][X][d1][d2]... with permute_order = "zyx", created on first write
template <typename T>
void Reader::WriteSlab(const T *data, const std::vector<int> &extra_dims, const char *dset_name)
{
    if (results_file_id < 0)
        throw std::runtime_error("WriteSlab: results file is not open");
    const Eigen::Index   extra = std::accumulate(extra_dims.begin(), extra_dims.end(), Eigen::Index(1), std::multiplies<>());
    std::vector<hsize_t> count{hsize_t(dims[2]), hsize_t(dims[1]), hsize_t(local_n0)}, offset(3 + extra_dims.size(), 0);
    count.insert(count.end(), extra_dims.begin(), extra_dims.end());
    offset[2] = local_0_start;

    std::unique_ptr<T[]> zyx(new T[local_n0 * dims[1] * dims[2] * extra]);
    swap_xz(data, zyx.get(), local_n0, dims[1], dims[2], extra);

    safe_create_group(results_file_id, dset_name);
    hid_t dset_id;
    if (H5Lexists(results_file_id, dset_name, H5P_DEFAULT) > 0) {
        dset_id = H5Dopen2(results_file_id, dset_name, H5P_DEFAULT);
    } else {
        std::vector<hsize_t> shape(count);
        shape[2]    = dims[0]; // all of X
        hid_t space = H5Screate_simple(shape.size(), shape.data(), nullptr);
        dset_id     = H5Dcreate2(results_file_id, dset_name, h5_native_type<T>(), space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t str = H5Tcopy(H5T_C_S1), scalar = H5Screate(H5S_SCALAR);
        H5Tset_size(str, 4);
        hid_t attr = H5Acreate2(dset_id, "permute_order", str, scalar, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, str, "zyx");
        H5Aclose(attr);
        H5Sclose(scalar);
        H5Tclose(str);
        H5Sclose(space);
    }

    hid_t filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset.data(), nullptr, count.data(), nullptr);
    hid_t memspace = H5Screate_simple(count.size(), count.data(), nullptr);
    hid_t xfer     = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(xfer, H5FD_MPIO_COLLECTIVE);
    const herr_t status = H5Dwrite(dset_id, h5_native_type<T>(), memspace, filespace, xfer, zyx.get());
    H5Pclose(xfer);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dset_id);
    if (status < 0)
        throw std::runtime_error(string("WriteSlab: cannot write ") + dset_name);
}

template <typename T>
std::vector<T> Reader::ReadData(const string &file, const string &dset_name)
{
    hid_t          file_id;
    hid_t          dset_id = h5_open_dataset(file, dset_name, file_id);
    hid_t          space   = H5Dget_space(dset_id);
    std::vector<T> data(H5Sget_simple_extent_npoints(space));
    const herr_t   status = H5Dread(dset_id, h5_native_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    H5Sclose(space);
    H5Dclose(dset_id);
    H5Fclose(file_id);
    if (status < 0)
        throw std::runtime_error("ReadData: cannot read " + dset_name + " in " + file);
    return data;
}

// This rank's slab of a per-voxel dataset on the grid, [Z][Y][X][d1][d2]...
// (or [X][Y][Z]..., see h5_is_zyx), into data as [X][Y][Z][d1][d2]...
template <typename T>
void Reader::ReadSlab(T *data, const std::vector<int> &extra_dims, const string &file, const string &dset_name) const
{
    hid_t      file_id;
    hid_t      dset_id = h5_open_dataset(file, dset_name, file_id);
    const bool zyx     = h5_is_zyx(dset_id);

    const Eigen::Index   extra = std::accumulate(extra_dims.begin(), extra_dims.end(), Eigen::Index(1), std::multiplies<>());
    std::vector<hsize_t> count, offset; // all of Y and Z, this rank's part of X
    if (zyx) {
        count  = {hsize_t(dims[2]), hsize_t(dims[1]), hsize_t(local_n0)};
        offset = {0, 0, hsize_t(local_0_start)};
    } else {
        count  = {hsize_t(local_n0), hsize_t(dims[1]), hsize_t(dims[2])};
        offset = {hsize_t(local_0_start), 0, 0};
    }
    count.insert(count.end(), extra_dims.begin(), extra_dims.end());
    offset.resize(count.size(), 0);

    hid_t filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset.data(), nullptr, count.data(), nullptr);
    hid_t                memspace = H5Screate_simple(count.size(), count.data(), nullptr);
    std::unique_ptr<T[]> zyx_data(zyx ? new T[local_n0 * dims[1] * dims[2] * extra] : nullptr);
    const herr_t         status = H5Dread(dset_id, h5_native_type<T>(), memspace, filespace, H5P_DEFAULT, zyx ? zyx_data.get() : data);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dset_id);
    H5Fclose(file_id);
    if (status < 0)
        throw std::runtime_error("ReadSlab: cannot read " + dset_name + " in " + file + " on the microstructure grid");
    if (zyx)
        swap_xz(zyx_data.get(), data, dims[2], dims[1], local_n0, extra);
}

template <typename T>
void Reader::writeData(const char *fieldName, int load_idx, int time_idx, const T *data, const hsize_t *shape, int rank)
{
    if (std::find(resultsToWrite.begin(), resultsToWrite.end(), fieldName) == resultsToWrite.end()) {
        return;
    }
    char name[5096];
    snprintf(name, sizeof(name), "%s/load%i/time_step%i/%s", dataset_name, load_idx, time_idx, fieldName);
    WriteData(data, name, shape, rank);
}

template <typename T>
void Reader::writeSlab(const char *fieldName, int load_idx, int time_idx, const T *data, const std::vector<int> &extra_dims)
{
    if (std::find(resultsToWrite.begin(), resultsToWrite.end(), fieldName) == resultsToWrite.end()) {
        return;
    }
    char name[5096];
    snprintf(name, sizeof(name), "%s/load%i/time_step%i/%s", dataset_name, load_idx, time_idx, fieldName);
    WriteSlab(data, extra_dims, name);
}

#endif
