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

// Constructor
MicroSimulation::MicroSimulation(int sim_id){
    // If used with the Micro Manager, MPI cannot be initialized again but
    // if the python bindings are used standalone, MPI should be initialized
    #ifdef USE_MPI
    MPI_Init(NULL, NULL);
    #endif
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // initialize fftw mpi
    fftw_mpi_init();

    // Convert the input file path to char* and read the input file
    const char* input_files_path = "input.json";
    int input_files_path_length = strlen(input_files_path) + 1;
    in_place_temp_path = new char[input_files_path_length];
    strcpy(in_place_temp_path, input_files_path);
    reader.ReadInputFile(in_place_temp_path);

    reader.ReadMS(3);

    // 3 signifies mechanics problems (maybe?) check with Sanath
    matmodel = createMatmodel<3>(reader);
    solver   = createSolver(reader, matmodel);
}

// Solve
py::dict MicroSimulation::solve(py::dict macro_data, double dt)
{
    // Create a pybind style Numpy array from macro_write_data["micro_vector_data"], which is a Numpy array
    py::array_t<double> macro_vector_data = macro_data["strains"].cast<py::array_t<double>>();

    // convert numpy array to std::vector.
    std::vector<double> g0_all = std::vector<double>(macro_vector_data.data(), macro_vector_data.data() + macro_vector_data.size());

    uint n_loads = g0_all.size() / matmodel->n_str;

    if (g0_all.size() % matmodel->n_str != 0)
        throw invalid_argument("Invalid length of loading g0");

    vector<double> g0(matmodel->n_str);

    matmodel->setGradient(g0);
    solver->solve();
    homogenized_stress = solver->get_homogenized_stress();

    // Get the stress from the solver
    std::vector<ssize_t> size_stress(4);
    size_stress[0] = reader.dims[0];
    size_stress[1] = reader.dims[1];
    size_stress[2] = reader.dims[2];
    size_stress[3] = matmodel->n_str;

    // Convert the stress array to a py::array_t<double>
    py::buffer_info info_stress{
        homogenized_stress,
        sizeof(double),
        py::format_descriptor<double>::format(),
        4,
        size_stress,
        {sizeof(double) * size_stress[1] * size_stress[2] * size_stress[3],
         sizeof(double) * size_stress[2] * size_stress[3], sizeof(double) * size_stress[3], sizeof(double)}
    };
    py::array_t<double> stress_array(info_stress);

    // Numerically calculate the stiffness matrix


    // Convert data to a py::dict again to send it back to the Micro Manager
    py::dict micro_write_data;

    // Add data to micro_write_data
    micro_write_data["stress"] = stress_array;
    micro_write_data["effective_stress"] = average_stress;

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
