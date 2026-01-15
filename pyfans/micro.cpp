// Micro simulation for mechanical problems
// In this file we solve a micro problem with FANS which is controlled by the Micro Manager
// This file is compiled with pybind11 to be available as a python module
//

#include "micro.hpp"
#include "setup.h"
#include "matmodel.h"
#include "mpi.h"

py::array_t<double> merge_arrays(py::array_t<double> array1, py::array_t<double> array2)
{
    // Ensure arrays are contiguous for efficient merging
    array1 = array1.attr("copy")();
    array2 = array2.attr("copy")();

    // Get numpy concatenate function
    py::object np          = py::module::import("numpy");
    py::object concatenate = np.attr("concatenate");

    // Concatenate the two arrays
    py::tuple           arrays = py::make_tuple(array1, array2);
    py::array_t<double> result = concatenate(arrays, py::int_(0)).cast<py::array_t<double>>();

    return result;
}

MicroSimulation::MicroSimulation(int sim_id, bool late_init, char *input_file)
    : _sim_id(sim_id)
{
    // initialize fftw mpi
    fftw_mpi_init();

    if (not late_init || sim_id >= 0) {
        // Input file name is hardcoded. TODO: Make it configurable
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

MicroSimulation::~MicroSimulation()
{
    Log::finalize();
}

py::dict MicroSimulation::solve(py::dict macro_data, double dt)
{
    const bool is_small_strain = std::holds_alternative<MaterialManager<3, 6> *>(matmanager);
    // Time step value dt is not used currently, but is available for future use

    // Create a pybind style Numpy array from macro_write_data["micro_vector_data"], which is a Numpy array
    py::array_t<double> strain1 = macro_data["strains1to3"].cast<py::array_t<double>>();
    py::array_t<double> strain2 = macro_data["strains4to6"].cast<py::array_t<double>>();

    py::array_t<double> strain = merge_arrays(strain1, strain2);
    if (not is_small_strain) {
        py::array_t<double> strain3 = macro_data["strains7to9"].cast<py::array_t<double>>();
        strain                      = merge_arrays(strain, strain3);
    }

    std::vector<double> g0 = std::vector<double>(strain.data(), strain.data() + strain.size()); // convert numpy array to std::vector.

    VectorXd homogenized_stress;

    std::visit([&](auto &mm) { mm->set_gradient(g0); }, matmanager);

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
        micro_write_data["stresses1to3"] = make_py_array(homogenized_stress[0], homogenized_stress[1], homogenized_stress[2]);
        micro_write_data["stresses4to6"] = make_py_array(homogenized_stress[3], homogenized_stress[4], homogenized_stress[5]);
        micro_write_data["cmat1"]        = make_py_array(C(0, 0), C(0, 1), C(0, 2));
        micro_write_data["cmat2"]        = make_py_array(C(0, 3), C(0, 4), C(0, 5));
        micro_write_data["cmat3"]        = make_py_array(C(1, 1), C(1, 2), C(1, 3));
        micro_write_data["cmat4"]        = make_py_array(C(1, 4), C(1, 5), C(2, 2));
        micro_write_data["cmat5"]        = make_py_array(C(2, 3), C(2, 4), C(2, 5));
        micro_write_data["cmat6"]        = make_py_array(C(3, 3), C(3, 4), C(3, 5));
        micro_write_data["cmat7"]        = make_py_array(C(4, 4), C(4, 5), C(5, 5));
    } else {
        micro_write_data["stresses1to3"] = make_py_array(homogenized_stress[0], homogenized_stress[1], homogenized_stress[2]);
        micro_write_data["stresses4to6"] = make_py_array(homogenized_stress[3], homogenized_stress[4], homogenized_stress[5]);
        micro_write_data["stresses7to9"] = make_py_array(homogenized_stress[6], homogenized_stress[7], homogenized_stress[8]);
    }

    return micro_write_data;
}

py::array_t<int> to_int_array(const std::vector<int> &v)
{
    py::array_t<int> arr(v.size());
    std::memcpy(arr.mutable_data(), v.data(), v.size() * sizeof(int));
    return arr;
}

py::dict MicroSimulation::get_state()
{
    py::dict state;

    Log::io->debug() << "Serializing reader\n";
    Serializable::buffer_t reader_state = reader.serialize_full();
    state["reader_state"]               = to_int_array(reader_state.as_int_vector());
    // Log::io->debug() << "Serializing solver\n";
    // Serializable::buffer_t solver_state = std::visit([](auto &s) { return s->serialize_full(); }, solver);
    // state["solver_state"] = to_int_array(solver_state.as_int_vector());
    // Log::io->debug() << "Serializing material manager\n";
    // Serializable::buffer_t matmngr_state = std::visit([](auto &m) { return m->serialize_full(); }, matmanager);
    // state["matmanager_state"] = to_int_array(matmngr_state.as_int_vector());
    Log::io->debug() << "Serializing micro sim done\n";

    return state;
}

void MicroSimulation::set_state(py::dict &state)
{
    auto                   py_reader_state = state["reader_state"].cast<py::array_t<int>>();
    std::vector<int>       reader_state_i(py_reader_state.data(), py_reader_state.data() + py_reader_state.size());
    Serializable::buffer_t reader_state = Serializable::buffer_t::from_int_vector(std::move(reader_state_i));

    // auto py_solver_state = state["solver_state"].cast<py::array_t<int>>();
    // std::vector<int> solver_state_i(py_solver_state.data(), py_solver_state.data() + py_solver_state.size());
    // Serializable::buffer_t solver_state = Serializable::buffer_t::from_int_vector(std::move(solver_state_i));

    // auto py_matmngr_state = state["matmanager_state"].cast<py::array_t<int>>();
    // std::vector<int> matmngr_state_i(py_matmngr_state.data(), py_matmngr_state.data() + py_matmngr_state.size());
    // Serializable::buffer_t matmngr_state = Serializable::buffer_t::from_int_vector(std::move(matmngr_state_i));

    reader.deserialize_full(reader_state);
    if (reader.strain_type == "small") {
        auto *mat_ptr = createMaterialManager<3, 6>(reader);
        // mat_ptr->deserialize_full(matmngr_state);
        matmanager    = mat_ptr;
        auto *sol_ptr = createSolver<3, 6>(reader, mat_ptr);
        // sol_ptr->init_fundamentalSolutionBuffer();
        // sol_ptr->matmanager = std::get<MaterialManager<3, 6> *>(matmanager);
        // sol_ptr->deserialize_full(solver_state);
        solver = sol_ptr;

    } else {
        auto *mat_ptr = createMaterialManager<3, 9>(reader);
        // mat_ptr->deserialize_full(matmngr_state);
        matmanager    = mat_ptr;
        auto *sol_ptr = createSolver<3, 9>(reader, mat_ptr);
        // sol_ptr->init_fundamentalSolutionBuffer();
        // sol_ptr->matmanager = std::get<MaterialManager<3, 9> *>(matmanager);
        // sol_ptr->deserialize_full(solver_state);
        solver = sol_ptr;
    }
}

int MicroSimulation::get_id()
{
    return _sim_id;
}

PYBIND11_MODULE(PyFANS, m)
{
    // optional docstring
    m.doc() = "FANS for Micro Manager";

    py::class_<MicroSimulation>(m, "MicroSimulation")
        .def(py::init<int>())
        .def("solve", &MicroSimulation::solve)
        .def(py::pickle(
            [](MicroSimulation &m) {
                return py::make_tuple(m.get_id(), m.get_state());
            },
            [](py::tuple t) {
                int      id    = t[0].cast<int>();
                py::dict state = t[1].cast<py::dict>();

                auto m = std::make_unique<MicroSimulation>(id, true);
                m->set_state(state);
                return m;
            }))
        .def("set_state", &MicroSimulation::set_state)
        .def("get_state", &MicroSimulation::get_state)
        .def("get_global_id", &MicroSimulation::get_id);
}
