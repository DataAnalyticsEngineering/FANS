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

 // Input file path has to be given in file path
  const char* input_path = "/home/torben/Documents/git/FANS/test/input_files/sphere_ThermalLinear.json";
  int input_path_length = strlen(input_path) + 1;  // Add 1 for null terminator

  char* input_temp_path = new char[input_path_length];  // Allocate memory for the path
  strcpy(input_temp_path, input_path);   

    // Input file path has to be given in file path
  const char* output_path = "/home/torben/Documents/git/FANS/test/input_files/sphere_ThermalLinear.json";
  int out_path_length = strlen(output_path) + 1;  // Add 1 for null terminator

  out_temp_path = new char[out_path_length];  // Allocate memory for the path
  strcpy(out_temp_path, output_path);   


    Reader reader;
    reader.ReadInputFile(input_temp_path);

    if(reader.problemType == "thermal"){

        for (int i = 0; i < world_size; i++){
    	    if(i == world_rank){
            	reader.ReadMS(1);
      	    }
      	    MPI_Barrier(MPI_COMM_WORLD);
        }
        reader.ComputeVolumeFractions();

        

        if(reader.matmodel == "ThermalLinear"){
            matmodel = new ThermalLinear(reader.l_e, reader.materialProperties);
        }else{
            throw invalid_argument(reader.matmodel + " is not a valid matmodel");
        }

        if(reader.method == "fp"){
            solver = new SolverFP<1>(reader, matmodel);
        }else if(reader.method == "cg"){
            solver = new SolverCG<1>(reader, matmodel);
        }else{
            throw invalid_argument(reader.method + " is not a valid method");
        }
    } else if(reader.problemType == "mechanical"){

        for (int i = 0; i < world_size; i++){
    	    if(i == world_rank){
            	reader.ReadMS(3);
      	    }
      	    MPI_Barrier(MPI_COMM_WORLD);
        }
        reader.ComputeVolumeFractions();


        if(reader.matmodel == "MechLinear"){
            matmodel = new MechLinear(reader.l_e, reader.materialProperties);
        }else if(reader.matmodel == "HyperElastic"){
            matmodel = new HyperElastic(reader.l_e, reader.materialProperties);
        }else{
            throw invalid_argument(reader.matmodel + " is not a valid matmodel");
        }

        if(reader.method == "fp"){
            solver = new SolverFP<3>(reader, matmodel);
        }else if(reader.method == "cg"){
            solver = new SolverCG<3>(reader, matmodel);
        }else{
            throw invalid_argument(reader.method + " is not a valid method");
        }
    }

    // add micro_scalar_data and micro_vector_data to micro_write_data
    micro_write_data["effective_stress"] = 0;
}

// Solve
py::dict MicroSimulation::solve(py::dict macro_data, double dt)
{

    // Create a pybind style Numpy array from macro_write_data["micro_vector_data"], which is a Numpy array
    py::array_t<double> macro_vector_data = macro_data["g0"].cast<py::array_t<double>>();
    _g0 = std::vector<double>(macro_vector_data.data(), macro_vector_data.data() + macro_vector_data.size()); // convert numpy array to std::vector.

    vector<double> g0_all = _g0;
    uint n_loads = g0_all.size() / matmodel->n_str;
    if (g0_all.size() % matmodel->n_str != 0) throw invalid_argument("Invalid length of loading g0");a

    vector<double> g0(matmodel->n_str);
    for(int i_load = 0; i_load < n_loads; i_load++){
        for (int i = 0; i < matmodel->n_str; ++i) {
            g0[i] = g0_all[i_load*matmodel->n_str + i];
        }

        matmodel->setGradient(g0);
        solver->solve();
        std::vector<double> average_stress = solver->postprocess(reader, out_temp_path, i_load);
        
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
    _micro_scalar_data = state[0].cast<double>();
    _state = state[1].cast<double>();
}

// This function needs to return variables which can fully define the state of a micro simulation
py::list MicroSimulation::get_state() const
{
    std::vector<double> state{_micro_scalar_data, _state};
    py::list state_python = py::cast(state);
    return state_python;
}

PYBIND11_MODULE(micro_dummy, m) {
    // optional docstring
    m.doc() = "pybind11 micro dummy plugin";

    py::class_<MicroSimulation>(m, "MicroSimulation")
        .def(py::init<int>())
        .def("initialize", &MicroSimulation::initialize)
        .def("solve", &MicroSimulation::solve)
        .def("get_state", &MicroSimulation::get_state)
        .def("set_state", &MicroSimulation::set_state)
        // Pickling support does not work currently, as there is no way to pass the simulation ID to the new instance ms.
        .def(py::pickle( // https://pybind11.readthedocs.io/en/latest/advanced/classes.html#pickling-support
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
            }
        ));
}
