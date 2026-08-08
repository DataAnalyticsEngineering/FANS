// Micro simulation for mechanical problems
// In this file we solve a micro problem with FANS which is controlled by the Micro Manager
// This file is compiled with pybind11 to be available as a python module
//

#include "micro.hpp"
#include "setup.h"
#include "matmodel.h"

PyFANSConfig load_config(const std::string &file_path)
{
    PyFANSConfig       config{};
    const PyFANSConfig defaults = config;

    try {
        ifstream i(file_path);
        json     j;
        i >> j;

        if (j.contains("no_mpi") && j["no_mpi"].get<bool>())
            config.disable_mpi = true;
        if (j.contains("log-level"))
            config.logging.level = spdlog::level::from_str(j.at("log-level").get<std::string>());
    } catch (const std::exception &e) {
        fprintf(stderr, "ERROR trying to read config file '%s' for pyFANS: %s\n", file_path.c_str(), e.what());
        fprintf(stderr, "Falling back to default values.\n");
        return defaults;
    }

    return config;
}

MicroSimulation::MicroSimulation(int sim_id, bool late_init, const std::string &input_file, const std::string &config_file)
    : _sim_id(sim_id), _input_file(input_file)
{
    // initialize fftw mpi
    fftw_mpi_init();

    const auto config = load_config(config_file);

    // Avoiding reader re-initialization due to unnecessary buffer copies
    reader.communicator = MPI_COMM_WORLD;
    if (config.disable_mpi)
        reader.communicator = MPI_COMM_SELF;
    MPI_Comm_rank(reader.communicator, &reader.world_rank);
    MPI_Comm_size(reader.communicator, &reader.world_size);

    Log::init(config.logging);

    if (not late_init || sim_id >= 0) {
        reader.ReadInputFile(input_file);
        reader.ReadMS(3);

        if (reader.strain_type == "small") {
            matmanager = createMaterialManager<3, 6>(reader);
            solver     = createSolver<3, 6>(reader, std::get<MaterialManager<3, 6> *>(matmanager));
        } else {
            matmanager = createMaterialManager<3, 9>(reader);
            solver     = createSolver<3, 9>(reader, std::get<MaterialManager<3, 9> *>(matmanager));
        }
    }
}

std::vector<double> merge_arrays(const std::vector<double> &v1, const std::vector<double> &v2)
{
    std::vector<double> res;
    res.resize(v1.size() + v2.size());
    std::copy(v1.begin(), v1.end(), res.begin());
    std::copy(v2.begin(), v2.end(), res.begin() + v1.size());
    return res;
}

std::vector<double> conv_to_vector(const py::array_t<double> &arr, const int size)
{
    std::vector<double> res;
    res.resize(size);
    // arr not necessarily contiguous
    for (int i = 0; i < size; i++)
        res[i] = arr.at(i);
    return res;
}

py::dict MicroSimulation::solve(const py::dict &macro_data, double dt)
{
    const bool is_small_strain = std::holds_alternative<MaterialManager<3, 6> *>(matmanager);
    // Time step value dt is not used currently, but is available for future use

    // Create a pybind style Numpy array from macro_write_data["micro_vector_data"], which is a Numpy array
    std::vector<double> strain1 = conv_to_vector(macro_data["Strains1to3"].cast<py::array_t<double>>(), 3);
    std::vector<double> strain2 = conv_to_vector(macro_data["Strains4to6"].cast<py::array_t<double>>(), 3);

    std::vector<double> strain = merge_arrays(strain1, strain2);
    if (not is_small_strain) {
        std::vector<double> strain3 = conv_to_vector(macro_data["Strains7to9"].cast<py::array_t<double>>(), 3);
        strain                      = merge_arrays(strain, strain3);
    }

    VectorXd homogenized_stress;

    std::visit([&](auto &mm) { mm->set_gradient(strain); }, matmanager);

    std::visit([](auto &s) { s->solve(); }, solver);

    homogenized_stress = std::visit([](auto &s) -> VectorXd { return s->get_homogenized_stress(); }, solver);

    MatrixXd C = std::visit([&](auto &s) -> MatrixXd { return s->get_homogenized_tangent(pert_param); }, solver);

    // Convert data to a py::dict again to send it back to the Micro Manager
    py::dict micro_write_data;

    auto make_py_array = [](double d0, double d1, double d2) {
        py::array_t<double> arr(3);
        auto                buf = arr.mutable_unchecked<1>();
        buf(0)                  = d0;
        buf(1)                  = d1;
        buf(2)                  = d2;
        return arr;
    };

    // Add stress and stiffness matrix data to Python dict to be returned
    if (is_small_strain) {
        micro_write_data["Stresses1to3"] = make_py_array(homogenized_stress[0], homogenized_stress[1], homogenized_stress[2]);
        micro_write_data["Stresses4to6"] = make_py_array(homogenized_stress[3], homogenized_stress[4], homogenized_stress[5]);
        micro_write_data["Cmat1"]        = make_py_array(C(0, 0), C(0, 1), C(0, 2));
        micro_write_data["Cmat2"]        = make_py_array(C(0, 3), C(0, 4), C(0, 5));
        micro_write_data["Cmat3"]        = make_py_array(C(1, 1), C(1, 2), C(1, 3));
        micro_write_data["Cmat4"]        = make_py_array(C(1, 4), C(1, 5), C(2, 2));
        micro_write_data["Cmat5"]        = make_py_array(C(2, 3), C(2, 4), C(2, 5));
        micro_write_data["Cmat6"]        = make_py_array(C(3, 3), C(3, 4), C(3, 5));
        micro_write_data["Cmat7"]        = make_py_array(C(4, 4), C(4, 5), C(5, 5));
    } else {
        micro_write_data["Stresses1to3"] = make_py_array(homogenized_stress[0], homogenized_stress[1], homogenized_stress[2]);
        micro_write_data["Stresses4to6"] = make_py_array(homogenized_stress[3], homogenized_stress[4], homogenized_stress[5]);
        micro_write_data["Stresses7to9"] = make_py_array(homogenized_stress[6], homogenized_stress[7], homogenized_stress[8]);
        micro_write_data["Cmat1"]        = make_py_array(C(0, 0), C(0, 1), C(0, 2));
        micro_write_data["Cmat2"]        = make_py_array(C(0, 3), C(0, 4), C(0, 5));
        micro_write_data["Cmat3"]        = make_py_array(C(0, 6), C(0, 7), C(0, 8));
        micro_write_data["Cmat4"]        = make_py_array(C(1, 0), C(1, 1), C(1, 2));
        micro_write_data["Cmat5"]        = make_py_array(C(1, 3), C(1, 4), C(1, 5));
        micro_write_data["Cmat6"]        = make_py_array(C(1, 6), C(1, 7), C(1, 8));
        micro_write_data["Cmat7"]        = make_py_array(C(2, 0), C(2, 1), C(2, 2));
        micro_write_data["Cmat8"]        = make_py_array(C(2, 3), C(2, 4), C(2, 5));
        micro_write_data["Cmat9"]        = make_py_array(C(2, 6), C(2, 7), C(2, 8));
        micro_write_data["Cmat10"]       = make_py_array(C(3, 0), C(3, 1), C(3, 2));
        micro_write_data["Cmat11"]       = make_py_array(C(3, 3), C(3, 4), C(3, 5));
        micro_write_data["Cmat12"]       = make_py_array(C(3, 6), C(3, 7), C(3, 8));
        micro_write_data["Cmat13"]       = make_py_array(C(4, 0), C(4, 1), C(4, 2));
        micro_write_data["Cmat14"]       = make_py_array(C(4, 3), C(4, 4), C(4, 5));
        micro_write_data["Cmat15"]       = make_py_array(C(4, 6), C(4, 7), C(4, 8));
        micro_write_data["Cmat16"]       = make_py_array(C(5, 0), C(5, 1), C(5, 2));
        micro_write_data["Cmat17"]       = make_py_array(C(5, 3), C(5, 4), C(5, 5));
        micro_write_data["Cmat18"]       = make_py_array(C(5, 6), C(5, 7), C(5, 8));
        micro_write_data["Cmat19"]       = make_py_array(C(6, 0), C(6, 1), C(6, 2));
        micro_write_data["Cmat20"]       = make_py_array(C(6, 3), C(6, 4), C(6, 5));
        micro_write_data["Cmat21"]       = make_py_array(C(6, 6), C(6, 7), C(6, 8));
        micro_write_data["Cmat22"]       = make_py_array(C(7, 0), C(7, 1), C(7, 2));
        micro_write_data["Cmat23"]       = make_py_array(C(7, 3), C(7, 4), C(7, 5));
        micro_write_data["Cmat24"]       = make_py_array(C(7, 6), C(7, 7), C(7, 8));
        micro_write_data["Cmat25"]       = make_py_array(C(8, 0), C(8, 1), C(8, 2));
        micro_write_data["Cmat26"]       = make_py_array(C(8, 3), C(8, 4), C(8, 5));
        micro_write_data["Cmat27"]       = make_py_array(C(8, 6), C(8, 7), C(8, 8));
    }

    return micro_write_data;
}

py::dict MicroSimulation::get_state()
{
    // TODO populate state
    py::dict state;
    return state;
}

void MicroSimulation::set_state(const py::dict &state)
{
    // TODO read from state, not file
    reader.FreeMS();
    reader.ReadInputFile(_input_file);
    reader.ReadMS(3);
    if (reader.strain_type == "small") {
        delete std::get<MaterialManager<3, 6> *>(matmanager);
        auto *mat_ptr = createMaterialManager<3, 6>(reader);
        matmanager    = mat_ptr;
        delete std::get<Solver<3, 6> *>(solver);
        auto *sol_ptr = createSolver<3, 6>(reader, mat_ptr);
        solver        = sol_ptr;

    } else {
        delete std::get<MaterialManager<3, 9> *>(matmanager);
        auto *mat_ptr = createMaterialManager<3, 9>(reader);
        matmanager    = mat_ptr;
        delete std::get<Solver<3, 9> *>(solver);
        auto *sol_ptr = createSolver<3, 9>(reader, mat_ptr);
        solver        = sol_ptr;
    }
}

int MicroSimulation::get_id()
{
    return _sim_id;
}

void MicroSimulation::set_id(const int id)
{
    _sim_id = id;
}

PYBIND11_MODULE(PyFANS, m)
{
    // optional docstring
    m.doc() = "FANS for Micro Manager";

    py::class_<MicroSimulation>(m, "MicroSimulation")
        .def(py::init<int, bool, const std::string &, const std::string &>(),
             py::arg("sim_id"),
             py::arg("late_init")  = false,
             py::arg("input_file") = "input.json",
             py::arg("config_file") = "pyfans-config.json")
        .def("solve", &MicroSimulation::solve, py::return_value_policy::automatic)
        .def("set_state", &MicroSimulation::set_state)
        .def("get_state", &MicroSimulation::get_state, py::return_value_policy::automatic)
        .def("get_global_id", &MicroSimulation::get_id)
        .def("set_global_id", &MicroSimulation::set_id);
}
