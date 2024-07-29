// Micro simulation
// In this file we solve a micro problem with FANS which is controlled by the Micro Manager
// This file is compiled with pybind11 to be available as a python module
//
// To check if python is able to import it, run:
// python3 -c "import micro; micro.MicroSimulation(1)"
// from the same directory

#include "micro.hpp"

// Constructor
MicroSimulation::MicroSimulation(int sim_id){
    #ifdef USE_MPI
    MPI_Init(NULL, NULL);
    #endif
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // initialize fftw mpi
    fftw_mpi_init();
    std::vector<double> average_stress;
    
    const char* input_path = "input.json";
    int input_path_length = strlen(input_path) + 1;  // Add 1 for null terminator

    input_temp_path = new char[input_path_length];  // Allocate memory for the path
    strcpy(input_temp_path, input_path);
    string input = reader.ReadFileLocations(input_temp_path);
    const char* input_files_path = input.c_str();
    int input_files_path_length = strlen(input_files_path) + 1;  // Add 1 for null terminator
    in_place_temp_path = new char[input_files_path_length];  // Allocate memory for the path
    strcpy(in_place_temp_path, input_files_path);   
    reader.ReadInputFile(in_place_temp_path);


    // from main.cpp
    if (reader.problemType == "thermal")
    {
        throw invalid_argument("Use the thermal simulation instead");
    }
    else if (reader.problemType == "mechanical")
    {

        for (int i = 0; i < world_size; i++)
        {
            if (i == world_rank)
            {
                reader.ReadMS(3);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        reader.ComputeVolumeFractions();


        if (reader.matmodel == "MechLinear")
        {
            matmodel = new MechLinear(reader.l_e, reader.materialProperties);
        }
        else if (reader.matmodel == "HyperElastic")
        {
            matmodel = new HyperElastic(reader.l_e, reader.materialProperties);
        }
        else
        {
            throw invalid_argument(reader.matmodel + " is not a valid matmodel");
        }

        if (reader.method == "fp")
        {
            solver = new SolverFP<3>(reader, matmodel);
        }
        else if (reader.method == "cg")
        {
            solver = new SolverCG<3>(reader, matmodel);
        }
        else
        {
            throw invalid_argument(reader.method + " is not a valid method");
        }
    }
}

// Solve
py::dict MicroSimulation::solve(py::dict macro_data, double dt)
{
    std::cout << "Solving micro problem" << std::endl;
    std::vector<double> average_stress;
    std::vector<double> average_strain;
    const char* output_path = reader.output_path.c_str();
    int out_path_length = strlen(output_path) + 1;  // Add 1 for null terminator

    out_temp_path = new char[out_path_length];  // Allocate memory for the path
    strcpy(out_temp_path, output_path); 
    // Create a pybind style Numpy array from macro_write_data["micro_vector_data"], which is a Numpy array
    py::array_t<double> macro_vector_data = macro_data["g0"].cast<py::array_t<double>>();
    std::vector<double> _g0 = std::vector<double>(macro_vector_data.data(), macro_vector_data.data() + macro_vector_data.size()); // convert numpy array to std::vector.

    vector<double> g0_all = _g0;

        uint n_loads = g0_all.size() / matmodel->n_str;

        if (g0_all.size() % matmodel->n_str != 0)
            throw invalid_argument("Invalid length of loading g0");
        vector<double> g0(matmodel->n_str);
        for (int i_load = 0; i_load < n_loads; i_load++)
        {
            for (int i = 0; i < matmodel->n_str; ++i)
            {
                g0[i] = g0_all[i_load * matmodel->n_str + i];
            }
            matmodel->setGradient(g0);
            solver->solve();
            std::tie(average_stress, average_strain) = solver->postprocess(reader, out_temp_path, i_load);
        }

    // Retrieve all values from postprocess and store them in a py::dict
    double *displacement = solver->v_u;
    std::vector<ssize_t> size_displacement(4);
    size_displacement[0] = reader.dims[0];
    size_displacement[1] = reader.dims[1];
    size_displacement[2] = reader.dims[2];
    size_displacement[3] = 3;

    // Convert the displacement and residual to a py::array_t<double>
    py::buffer_info disp_info{
        displacement, // Pointer to data (nullptr if the array is empty)
        sizeof(double), // Size of one scalar
        py::format_descriptor<double>::format(), // Python struct-style format descriptor
        4, // Number of dimensions
        size_displacement, // Buffer dimensions
        {sizeof(double) * size_displacement[1] * size_displacement[2] * size_displacement[3], sizeof(double) * size_displacement[2] * size_displacement[3], sizeof(double) * size_displacement[3], sizeof(double)} // Strides (in bytes) for each index
    };
    py::array_t<double> displacement_array(disp_info);

    double* stress = solver->stress;
    std::vector<ssize_t> size_stress(4);
    size_stress[0] = reader.dims[0];
    size_stress[1] = reader.dims[1];
    size_stress[2] = reader.dims[2];
    size_stress[3] = matmodel->n_str;

    py::buffer_info info_stress{
        stress, // Pointer to data (nullptr if the array is empty)
        sizeof(double), // Size of one scalar
        py::format_descriptor<double>::format(), // Python struct-style format descriptor
        4, // Number of dimensions
        size_stress, // Buffer dimensions
        {sizeof(double) * size_stress[1] * size_stress[2] * size_stress[3], sizeof(double) * size_stress[2] * size_stress[3], sizeof(double) * size_stress[3], sizeof(double)} // Strides (in bytes) for each index
    };
    py::array_t<double> stress_array(info_stress);

    // Convert data to a py::dict again to send it back to the Micro Manager
    py::dict micro_write_data;

    // add micro_scalar_data and micro_vector_data to micro_write_data

    micro_write_data["stress"] = stress_array;
    micro_write_data["effective_stress"] = average_stress;
    micro_write_data["displacement"] = displacement_array;
    return micro_write_data;
}


PYBIND11_MODULE(PyFANS, m)
{
    // optional docstring
    m.doc() = "pybind11 micro dummy plugin";

    py::class_<MicroSimulation>(m, "MicroSimulation")
        .def(py::init<int>())
        .def("solve", &MicroSimulation::solve);
}
