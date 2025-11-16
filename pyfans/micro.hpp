// This is the header file for the micro simulation class.
// It is included in the micro_cpp_dummy.cpp file and the micro_cpp_dummy.cpp file is compiled with pybind11 to create a python module.
// The python module is then imported in the Micro Manager.

#pragma once
#include <iostream>
#include <vector>

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h" // numpy arrays
#include "pybind11/stl.h"   // std::vector conversion

#include "general.h"
#include "matmodel.h"
#include "MaterialManager.h"
#include "solver.h"

namespace py = pybind11;

class MicroSimulation {
  public:
    MicroSimulation(int sim_id, char *input_file = "input.json");
    py::dict solve(py::dict macro_write_data, double dt);

  private:
    int    _sim_id;
    Reader reader;
    // Hardcoding mechanical models because these definitions need information from the input file.
    MaterialManager<3, 6> *matmanager;
    Solver<3, 6>          *solver;
    double                 pert_param = 1e-6; // scalar strain perturbation parameter
};
