// This is the header file for the micro simulation class.
// It is included in the micro_cpp_dummy.cpp file and the micro_cpp_dummy.cpp file is compiled with pybind11 to create a python module.
// The python module is then imported in the Micro Manager.

#pragma once
#include <iostream>
#include <vector>
#include <variant>

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h" // numpy arrays
#include "pybind11/stl.h"   // std::vector conversion

#include "general.h"
#include "matmodel.h"
#include "MaterialManager.h"
#include "solver.h"

namespace py = pybind11;

#ifndef PyFANS_CLASS_NAME
    #define PyFANS_CLASS_NAME MicroSimulation
#endif
#ifndef PyFANS_CLASS_NAME_STR
    #define PyFANS_CLASS_NAME_STR "MicroSimulation"
#endif
#ifndef PyFANS_INPUT_NAME
    #define PyFANS_INPUT_NAME "input.json"
#endif
#ifndef PyFANS_MODULE_NAME
    #define PyFANS_MODULE_NAME PyFANS
#endif


class PyFANS_CLASS_NAME {
  public:
    PyFANS_CLASS_NAME(int sim_id, bool late_init = false, char *input_file = PyFANS_INPUT_NAME);
    ~PyFANS_CLASS_NAME();
    py::dict solve(const py::dict &macro_write_data, double dt);

    py::dict get_state();
    void     set_state(const py::dict &state);

    int get_id();

  private:
    int    _sim_id;
    Reader reader;
    // Hardcoding mechanical models because these definitions need information from the input file.
    using matmanager_t = std::variant<MaterialManager<3, 6> *, MaterialManager<3, 9> *>;
    using solver_t     = std::variant<Solver<3, 6> *, Solver<3, 9> *>;
    matmanager_t matmanager;
    solver_t     solver;
    double       pert_param = 1e-6; // scalar strain perturbation parameter
};
