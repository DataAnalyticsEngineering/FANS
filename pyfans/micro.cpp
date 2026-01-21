// Micro simulation for mechanical problems
// In this file we solve a micro problem with FANS which is controlled by the Micro Manager
// This file is compiled with pybind11 to be available as a python module
//

#include "micro.hpp"
#include "setup.h"
#include "matmodel.h"
#include "mpi.h"

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
    // Log::finalize();
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
    for (int i = 0; i < size; i++) res[i] = arr.at(i);
    return res;
}

py::dict MicroSimulation::solve(const py::dict &macro_data, double dt)
{
    const bool is_small_strain = std::holds_alternative<MaterialManager<3, 6> *>(matmanager);
    // Time step value dt is not used currently, but is available for future use

    // Create a pybind style Numpy array from macro_write_data["micro_vector_data"], which is a Numpy array
    std::vector<double> strain1 = conv_to_vector(macro_data["strains1to3"].cast<py::array_t<double>>(), 3);
    std::vector<double> strain2 = conv_to_vector(macro_data["strains4to6"].cast<py::array_t<double>>(), 3);

    std::vector<double> strain = merge_arrays(strain1, strain2);
    if (not is_small_strain) {
        std::vector<double> strain3 = conv_to_vector(macro_data["strains7to9"].cast<py::array_t<double>>(), 3);
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

py::dict MicroSimulation::get_state()
{
    py::dict state;
    return state;
}

void MicroSimulation::set_state(const py::dict &state)
{
    reader.FreeMS();
    reader.ReadInputFile("input.json");
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

PYBIND11_MODULE(PyFANS, m)
{
    // optional docstring
    m.doc() = "FANS for Micro Manager";

    py::class_<MicroSimulation>(m, "MicroSimulation")
        .def(py::init<int>())
        .def("solve", &MicroSimulation::solve, py::return_value_policy::automatic)
        .def("set_state", &MicroSimulation::set_state)
        .def("get_state", &MicroSimulation::get_state, py::return_value_policy::automatic)
        .def("get_global_id", &MicroSimulation::get_id);
}
