// Micro simulation for mechanical problems
// In this file we solve a micro problem with FANS which is controlled by the Micro Manager
// This file is compiled with pybind11 to be available as a python module
//
// To check if python is able to import it, run:
// python3 -c "import micro; micro.MicroSimulation(1)"
// from the same directory

#include "micro.hpp"
#include "setup.h"
#include "matmodel.h"

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

// Constructor
MicroSimulation::MicroSimulation(int sim_id)
{
    MPI_Init(NULL, NULL);

    // initialize fftw mpi
    fftw_mpi_init();

    // Convert the input file path to char* and read the input file
    const char *input_files_path        = "input_files/test_LinearElastic.json";
    int         input_files_path_length = strlen(input_files_path) + 1;
    in_place_temp_path                  = new char[input_files_path_length];
    strcpy(in_place_temp_path, input_files_path);
    reader.ReadInputFile(in_place_temp_path);

    reader.ReadMS(3);

    // 3 signifies mechanics problems
    matmodel = createMatmodel<3>(reader);
    solver   = createSolver(reader, matmodel);
}

// Solve
py::dict MicroSimulation::solve(py::dict macro_data, double dt)
{
    std::cout << "Solving micro problem" << std::endl;

    // Create a pybind style Numpy array from macro_write_data["micro_vector_data"], which is a Numpy array
    py::array_t<double> strain1 = macro_data["strains1to3"].cast<py::array_t<double>>();
    py::array_t<double> strain2 = macro_data["strains4to6"].cast<py::array_t<double>>();

    py::array_t<double> strain = merge_arrays(strain1, strain2);
    std::vector<double> _g0    = std::vector<double>(strain.data(), strain.data() + strain.size()); // convert numpy array to std::vector.

    vector<double> g0_all = _g0;

    uint n_loads = g0_all.size() / matmodel->n_str;

    if (g0_all.size() % matmodel->n_str != 0)
        throw invalid_argument("Invalid length of loading g0");

    vector<double> g0(matmodel->n_str);

    VectorXd homogenized_stress;

    for (int i_load = 0; i_load < n_loads; i_load++) {
        for (int i = 0; i < matmodel->n_str; ++i) {
            g0[i] = g0_all[i_load * matmodel->n_str + i];
        }
        matmodel->setGradient(g0);
        solver->solve();
        homogenized_stress = solver->get_homogenized_stress();
    }

    auto C = solver->get_homogenized_tangent(pert_param);

    std::cout << "Homogenized stiffness matrix C: " << C << std::endl;

    // Convert data to a py::dict again to send it back to the Micro Manager
    py::dict micro_write_data;

    // Add stress and stiffness matrix data to Python dict to be returned
    std::vector<double> stress13     = {homogenized_stress[0], homogenized_stress[1], homogenized_stress[2]};
    micro_write_data["stresses1to3"] = stress13;
    std::vector<double> stress46     = {homogenized_stress[3], homogenized_stress[4], homogenized_stress[5]};
    micro_write_data["stresses4to6"] = stress46;
    std::vector<double> C_1          = {C(0, 0), C(0, 1), C(0, 2)};
    micro_write_data["cmat1"]        = C_1;
    std::vector<double> C_2          = {C(0, 3), C(0, 4), C(0, 5)};
    micro_write_data["cmat2"]        = C_2;
    std::vector<double> C_3          = {C(1, 1), C(1, 2), C(1, 3)};
    micro_write_data["cmat3"]        = C_3;
    std::vector<double> C_4          = {C(1, 4), C(1, 5), C(2, 2)};
    micro_write_data["cmat4"]        = C_4;
    std::vector<double> C_5          = {C(2, 3), C(2, 4), C(2, 5)};
    micro_write_data["cmat5"]        = C_5;
    std::vector<double> C_6          = {C(3, 3), C(3, 4), C(3, 5)};
    micro_write_data["cmat6"]        = C_6;
    std::vector<double> C_7          = {C(4, 4), C(4, 5), C(5, 5)};
    micro_write_data["cmat7"]        = C_7;

    return micro_write_data;
}

PYBIND11_MODULE(PyFANS, m)
{
    // optional docstring
    m.doc() = "FANS for Micro Manager";

    py::class_<MicroSimulation>(m, "MicroSimulation")
        .def(py::init<int>())
        .def("solve", &MicroSimulation::solve);
}
