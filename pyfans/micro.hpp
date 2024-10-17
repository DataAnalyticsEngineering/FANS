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
#include "solver.h"
#include "solverFP.h"
#include "solverCG.h"

namespace py = pybind11;

class MicroSimulation {
  public:
    MicroSimulation(int sim_id);
    py::dict solve(py::dict macro_write_data, double dt); // API according to the Micro Manager

  private:
    int                 _sim_id;
    std::vector<double> g0;
    double              _state;
    char               *input_temp_path;
    char               *in_place_temp_path;
    char               *out_temp_path;
    Reader              reader;
    Matmodel<3>        *matmodel;
    Solver<3>          *solver;

    double               pert_param = 1e-3; // scalar strain perturbation parameter
    std::vector<double>  homogenized_stress;
    std::vector<double>  unperturbed_stress;
    Matrix<double, 6, 6> C; // Homogenized stiffness matrix C
};