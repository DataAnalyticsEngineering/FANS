#include "general.h"
#include "reader.h"

#include <cstdlib>

Reader::Reader(const MPI_Comm &comm)
    : communicator(comm)
{
    MPI_Comm_rank(communicator, &world_rank);
    MPI_Comm_size(communicator, &world_size);
}

void Reader::ComputeVolumeFractions()
{
    const size_t local_size = local_n0 * dims[1] * dims[2];

    // The phases are 0 ... n_mat - 1
    int highest = *std::max_element(ms, ms + local_size);
    MPI_Allreduce(MPI_IN_PLACE, &highest, 1, MPI_INT, MPI_MAX, communicator);
    n_mat = highest + 1;

    // Voxels of each phase, summed over all ranks
    std::vector<long> voxels(n_mat, 0);
    for (size_t i = 0; i < local_size; i++)
        voxels[ms[i]]++;
    MPI_Allreduce(MPI_IN_PLACE, voxels.data(), n_mat, MPI_LONG, MPI_SUM, communicator);
    if (voxels[0] == 0)
        throw std::invalid_argument("Microstructure phase IDs must start at 0");

    Log::logger().info("# Number of materials: {} (from 0 to {})", n_mat, highest);
    Log::logger().info("# Volume fractions");
    for (int i = 0; i < n_mat; i++)
        Log::logger().info("# material {:4}    vol. frac. {:10.4f}%  ", i, 100. * voxels[i] / (double(dims[0]) * dims[1] * dims[2]));
}

void Reader ::ReadInputFile(const std::string &input_fn)
{
    try {
        ifstream i(input_fn);
        json     j;
        i >> j;
        inputJson = j; // Store complete input JSON for MaterialManager

        microstructure = j["microstructure"];
        ms_filename    = microstructure["filepath"].get<string>();
        ms_datasetname = microstructure["datasetname"].get<string>();
        if (ms_datasetname.empty())
            throw std::invalid_argument("datasetname must not be empty and must refer to a valid HDF5 path");
        if (ms_datasetname.front() != '/')
            ms_datasetname = "/" + ms_datasetname; // an absolute HDF5 path
        L = microstructure["L"].get<vector<double>>();

        // The group of the microstructure, and where the results go: "<ms_datasetname>_results/<results_prefix>"
        ms_group     = ms_datasetname.substr(0, ms_datasetname.find_last_of('/') + 1);
        dataset_name = ms_datasetname + "_results/" + j.value("results_prefix", "");

        errorParameters = j["error_parameters"];
        TOL             = errorParameters["tolerance"].get<double>();
        n_it            = j["n_it"].get<int>();

        extrapolate_displacement = j.value("extrapolate_displacement", extrapolate_displacement);

        if (j.contains("linesearch_parameters")) {
            ls_max_iter = j["linesearch_parameters"].value("max_iter", ls_max_iter);
            ls_tol      = j["linesearch_parameters"].value("tol", ls_tol);
            if (ls_max_iter < 1 || ls_tol <= 0.0)
                throw std::invalid_argument("linesearch_parameters: max_iter >= 1 and tol > 0 required");
        }

        problemType = j["problem_type"].get<string>();
        method      = j["method"].get<string>();

        strain_type = j.value("strain_type", "small");
        if (strain_type != "small" && strain_type != "large")
            throw std::invalid_argument("strain_type must be either 'small' or 'large'");

        FE_type = j.value("FE_type", "HEX8"); // full integration
        if (FE_type != "HEX8" && FE_type != "HEX8R" && FE_type != "BBAR")
            throw std::invalid_argument("FE_type must be one of: 'HEX8', 'HEX8R', or 'BBAR'");

        resultsToWrite = j["results"].get<vector<string>>(); // Read the results_to_write field

        load_cases.clear();
        const auto &ml = j["macroscale_loading"];
        if (!ml.is_array())
            throw std::runtime_error("macroscale_loading must be an array");

        // Determine the size of loading vector based on problem type and strain formulation
        int n_str;
        if (problemType == "thermal") {
            n_str = 3; // Temperature gradient components
        } else if (strain_type == "large") {
            n_str = 9; // Deformation gradient components (F11, F12, F13, F21, F22, F23, F31, F32, F33)
        } else {
            n_str = 6; // Small strain components (eps11, eps22, eps33, eps12, eps13, eps23)
        }

        for (const auto &entry : ml) {
            LoadCase lc;
            if (entry.is_array()) { // ---------- legacy pure-strain ----------
                lc.mixed   = false;
                lc.g0_path = entry.get<vector<vector<double>>>();
                lc.n_steps = lc.g0_path.size();
                if (lc.g0_path[0].size() != static_cast<size_t>(n_str))
                    throw std::invalid_argument("Invalid length of loading vector: expected " +
                                                std::to_string(n_str) + " components but got " +
                                                std::to_string(lc.g0_path[0].size()));
            } else { // ---------- mixed BC object ------------
                lc.mixed   = true;
                lc.mbc     = MixedBC::from_json(entry, n_str);
                lc.n_steps = lc.mbc.F_E_path.rows();
            }
            load_cases.push_back(std::move(lc));
        }

        // "time_step": one size for all steps, or a list per load case like macroscale_loading
        const json time_step = j.value("time_step", json(1.0));
        for (size_t c = 0; c < load_cases.size(); ++c) {
            auto &lc = load_cases[c];
            lc.dt    = time_step.is_number() ? vector<double>(lc.n_steps, time_step.get<double>()) : time_step.at(c).get<vector<double>>();
            if (lc.dt.size() != lc.n_steps || *std::min_element(lc.dt.begin(), lc.dt.end()) <= 0.0)
                throw std::invalid_argument("time_step of load case " + std::to_string(c + 1) + " must hold " +
                                            std::to_string(lc.n_steps) + " positive values");
        }

        Log::logger().info("# microstructure file name: \t '{}'", ms_filename);
        Log::logger().info("# microstructure dataset name: \t '{}'", ms_datasetname);
        Log::logger().info("# strain type: \t {}", strain_type);
        Log::logger().info("# problem type: \t {}", problemType);
        Log::logger().info("# FE type: \t {}", FE_type);
        Log::logger().info("# FANS error measure: \t {} {} error  ", errorParameters["type"].get<string>(), errorParameters["measure"].get<string>());
        Log::logger().info("# FANS Tolerance: \t {:10.5e}", errorParameters["tolerance"].get<double>());
        Log::logger().info("# Max iterations: \t {:6}", n_it);
    } catch (const std::exception &e) {
        Log::logger().error("ERROR trying to read input file '{}' for FANS: {}", input_fn, e.what());
        exit(10);
    }
}

// Creates the groups on the path of `name` that do not exist yet, e.g. /a and /a/b for /a/b/c
void Reader::safe_create_group(hid_t file, const string &path)
{
    for (size_t i = path.find('/', 1); i != string::npos; i = path.find('/', i + 1)) {
        const string group = path.substr(0, i);
        if (H5Lexists(file, group.c_str(), H5P_DEFAULT) <= 0)
            H5Gclose(H5Gcreate(file, group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }
}

void Reader ::ReadMS(int hm)
{
    // Grid size, and whether the file is laid out [Z][Y][X] (the default) or [X][Y][Z]
    hid_t   file_id;
    hid_t   dset_id = h5_open_dataset(ms_filename, ms_datasetname, file_id);
    hid_t   space   = H5Dget_space(dset_id);
    hsize_t file_dims[3];
    if (H5Sget_simple_extent_ndims(space) != 3)
        throw std::runtime_error("Microstructure dataset " + ms_datasetname + " must be 3D");
    H5Sget_simple_extent_dims(space, file_dims, nullptr);
    H5Sclose(space);

    is_zyx = h5_is_zyx(dset_id);
    H5Dclose(dset_id);
    H5Fclose(file_id);
    if (is_zyx) {
        Log::logger().info("# Using Z-Y-X dimension ordering for the microstructure data");
        dims = {int(file_dims[2]), int(file_dims[1]), int(file_dims[0])};
    } else {
        Log::logger().info("# Using X-Y-Z dimension ordering for the microstructure data");
        dims = {int(file_dims[0]), int(file_dims[1]), int(file_dims[2])};
    }
    l_e = {L[0] / dims[0], L[1] / dims[1], L[2] / dims[2]};

    Log::logger().info("# Grid size set to [{} x {} x {}] --> {} voxels", dims[0], dims[1], dims[2], dims[0] * dims[1] * dims[2]);
    Log::logger().info("# Microstructure length: [{:3.6f} x {:3.6f} x {:3.6f}]", L[0], L[1], L[2]);
    Log::logger().info("# Voxel length: [{:.8f}, {:.8f}, {:.8f}]", l_e[0], l_e[1], l_e[2]);

    if (dims[0] % 2 != 0)
        Log::logger().warn("[ FANS3D_Grid ] WARNING: n_x is not a multiple of 2");
    if (dims[1] % 2 != 0)
        Log::logger().warn("[ FANS3D_Grid ] WARNING: n_y is not a multiple of 2");
    if (dims[2] % 2 != 0)
        Log::logger().warn("[ FANS3D_Grid ] WARNING: n_z is not a multiple of 2");
    if (dims[0] / 4 < world_size)
        throw std::runtime_error("[ FANS3D_Grid ] ERROR: Please decrease the number of processes or increase the grid size to ensure that each process has at least 4 boxels in the x direction.");

    // This rank's part of the grid: x-slab [local_0_start, local_0_start + local_n0), see
    // https://www.fftw.org/fftw3_doc/Transposed-distributions.html
    const ptrdiff_t n[3] = {dims[0], dims[1], dims[2] / 2 + 1};
    alloc_local          = fftw_mpi_local_size_many_transposed(3, n, hm, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, communicator,
                                                               &local_n0, &local_0_start, &local_n1, &local_1_start);
    if (local_n0 < 4)
        throw std::runtime_error("[ FANS3D_Grid ] ERROR: Number of voxels in x-direction is less than 4 in process " + to_string(world_rank));

    ms = FANS_malloc<unsigned short>(size_t(local_n0) * dims[1] * dims[2]);
    ReadSlab(ms, {}, ms_filename, ms_datasetname);

    this->ComputeVolumeFractions();
}

hid_t h5_open_dataset(const string &file, const string &dset_name, hid_t &file_id)
{
    file_id = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0)
        throw std::runtime_error("Failed to open the HDF5 file " + file);
    hid_t dset_id = H5Dopen2(file_id, dset_name.c_str(), H5P_DEFAULT);
    if (dset_id < 0) {
        H5Fclose(file_id);
        throw std::runtime_error("No dataset " + dset_name + " in " + file);
    }
    return dset_id;
}

bool h5_is_zyx(hid_t dset_id)
{
    if (H5Aexists(dset_id, "permute_order") <= 0)
        return true;
    // a variable-length string (h5py's default) or a fixed-length one (C, MATLAB, WriteSlab)
    hid_t  attr_id   = H5Aopen(dset_id, "permute_order", H5P_DEFAULT);
    hid_t  attr_type = H5Aget_type(attr_id);
    string permute_order;
    if (H5Tis_variable_str(attr_type) > 0) {
        char *str = nullptr;
        if (H5Aread(attr_id, attr_type, &str) >= 0 && str != nullptr)
            permute_order = str;
        H5free_memory(str);
    } else {
        permute_order.resize(H5Tget_size(attr_type));
        if (H5Aread(attr_id, attr_type, permute_order.data()) < 0)
            permute_order.clear();
    }
    H5Tclose(attr_type);
    H5Aclose(attr_id);
    return permute_order.empty() || permute_order[0] == 'z' || permute_order[0] == 'Z';
}

std::vector<hsize_t> Reader::DataShape(const string &file, const string &dset_name)
{
    hid_t                file_id;
    hid_t                dset_id = h5_open_dataset(file, dset_name, file_id);
    hid_t                space   = H5Dget_space(dset_id);
    std::vector<hsize_t> shape(H5Sget_simple_extent_ndims(space));
    H5Sget_simple_extent_dims(space, shape.data(), nullptr);
    H5Sclose(space);
    H5Dclose(dset_id);
    H5Fclose(file_id);
    return shape;
}

void Reader::FreeMS()
{
    if (ms) {
        FANS_free(ms);
        ms = nullptr;
    }
}

void Reader::OpenResultsFile(const char *output_fn)
{
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, communicator, MPI_INFO_NULL);
    results_file_id = H5Fcreate(output_fn, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    if (results_file_id < 0) {
        throw std::runtime_error("Failed to create results file");
    }
}

void Reader::CloseResultsFile()
{
    if (results_file_id >= 0) {
        H5Fclose(results_file_id);
        results_file_id = -1;
    }
}

Reader::~Reader()
{
    FreeMS();
}
