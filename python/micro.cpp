// Micro simulation
// In this file we solve a micro problem with FANS which is controlled by the Micro Manager
// This dummy is written in C++ and is controllable via Python using pybind11
//
// Compile your pybind-11 wrapped code with:
//
// c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) micro.cpp -o micro$(python3-config --extension-suffix)
//
// To check if python is able to import it, run:
// python3 -c "import micro; micro.MicroSimulation(1)"
// from the same directory

#include "micro.hpp"

// Constructor
MicroSimulation::MicroSimulation(int sim_id) : _sim_id(sim_id), _state(0) {}

// Initialize
void MicroSimulation::initialize()
{

    MPI_Init(NULL, NULL);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // initialize fftw mpi
    fftw_mpi_init();

    // Set inputs manually for this test case (taken from sphere_ThermalLinear.json)
    std::strncpy(reader.ms_filename, "microstructures/sphere32.h5", sizeof(reader.ms_filename) - 1);
    std::strncpy(reader.ms_datasetname, "/sphere/32x32x32/ms", sizeof(reader.ms_datasetname) - 1);

    reader.L = {1.0, 1.0, 1.0};

    reader.problemType = "thermal";
    reader.matmodel = "ThermalLinear";

    reader.materialProperties["conductivity"] = {1, 100};

    reader.method = "cg";
    reader.TOL = 1e-6;
    reader.n_it = 100;

    reader.g0 = {
        0.0002, -0.0002, 0, 
        0.0004, -0.0004, 0, 
        0.0006, -0.0006, 0, 
        0.0008, -0.0008, 0, 
        0.001, -0.001, 0, 
        0.0012, -0.0012, 0, 
        0.0014, -0.0014, 0, 
        0.0016, -0.0016, 0, 
        0.0018, -0.0018, 0, 
        0.002, -0.002, 0, 
        0.0022, -0.0022, 0, 
        0.0024, -0.0024, 0, 
        0.0026, -0.0026, 0, 
        0.0028, -0.0028, 0, 
        0.003, -0.003, 0, 
        0.0032, -0.0032, 0, 
        0.0034, -0.0034, 0, 
        0.0036, -0.0036, 0, 
        0.0038, -0.0038, 0, 
        0.004, -0.004, 0, 
        0.0042, -0.0042, 0, 
        0.0044, -0.0044, 0, 
        0.0046, -0.0046, 0, 
        0.0048, -0.0048, 0, 
        0.005, -0.005, 0, 
        0.0052, -0.0052, 0, 
        0.0054, -0.0054, 0, 
        0.0056, -0.0056, 0, 
        0.0058, -0.0058, 0, 
        0.006, -0.006, 0, 
        0.0062, -0.0062, 0, 
        0.0064, -0.0064, 0, 
        0.0066, -0.0066, 0, 
        0.0068, -0.0068, 0, 
        0.007, -0.007, 0, 
        0.0072, -0.0072, 0, 
        0.0074, -0.0074, 0, 
        0.0076, -0.0076, 0, 
        0.0078, -0.0078, 0, 
        0.008, -0.008, 0, 
        0.0082, -0.0082, 0, 
        0.0084, -0.0084, 0, 
        0.0086, -0.0086, 0, 
        0.0088, -0.0088, 0, 
        0.009, -0.009, 0, 
        0.0092, -0.0092, 0, 
        0.0094, -0.0094, 0, 
        0.0096, -0.0096, 0, 
        0.0098, -0.0098, 0, 
        0.01, -0.01, 0
    };

    reader.resultsToWrite = {"microstructure", "stress", "strain", "displacement", "plastic_flag",
                             "stress_average", "strain_average", "absolute_error"};

    // from main.cpp
    if (reader.problemType == "thermal")
    {

        for (int i = 0; i < world_size; i++)
        {
            if (i == world_rank)
            {
                reader.ReadMS(1);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        reader.ComputeVolumeFractions();

        if (reader.matmodel == "ThermalLinear")
        {
            matmodel = new ThermalLinear(reader.l_e, reader.materialProperties);
        }
        else
        {
            throw invalid_argument(reader.matmodel + " is not a valid matmodel");
        }

        if (reader.method == "fp")
        {
            solver = new SolverFP<1>(reader, matmodel);
        }
        else if (reader.method == "cg")
        {
            solver = new SolverCG<1>(reader, matmodel);
        }
        else
        {
            throw invalid_argument(reader.method + " is not a valid method");
        }
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
            matmodel2 = new MechLinear(reader.l_e, reader.materialProperties);
        }
        else if (reader.matmodel == "HyperElastic")
        {
            matmodel2 = new HyperElastic(reader.l_e, reader.materialProperties);
        }
        else
        {
            throw invalid_argument(reader.matmodel + " is not a valid matmodel");
        }

        if (reader.method == "fp")
        {
            solver2 = new SolverFP<3>(reader, matmodel2);
        }
        else if (reader.method == "cg")
        {
            solver2 = new SolverCG<3>(reader, matmodel2);
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
    std::vector<double> average_stress;
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
        average_stress = solver->postprocess(reader, out_temp_path, i_load);
    }

    // Convert data to a py::dict again to send it back to the Micro Manager
    py::dict micro_write_data;

    // add micro_scalar_data and micro_vector_data to micro_write_data
    micro_write_data["effective_stress"] = average_stress;
    return micro_write_data;
}

// This function needs to set the complete state of a micro simulation
void MicroSimulation::set_state(py::list state)
{
}

// This function needs to return variables which can fully define the state of a micro simulation
py::list MicroSimulation::get_state() const
{
}

PYBIND11_MODULE(PyFANS, m)
{
    // optional docstring
    m.doc() = "pybind11 micro dummy plugin";

    py::class_<MicroSimulation>(m, "MicroSimulation")
        .def(py::init<int>())
        .def("initialize", &MicroSimulation::initialize)
        .def("solve", &MicroSimulation::solve)
        .def("get_state", &MicroSimulation::get_state)
        .def("set_state", &MicroSimulation::set_state)
        // Pickling support does not work currently, as there is no way to pass the simulation ID to the new instance ms.
        .def(py::pickle(                    // https://pybind11.readthedocs.io/en/latest/advanced/classes.html#pickling-support
            [](const MicroSimulation &ms) { // __getstate__
                return ms.get_state();
            },
            [](py::list t) { // __setstate__
                if (t.size() != 2)
                    throw std::runtime_error("Invalid state!");

                /* Create a new C++ instance */
                MicroSimulation ms(0);

                ms.set_state(t);

                return ms;
            }));
}
