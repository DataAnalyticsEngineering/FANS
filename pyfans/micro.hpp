// This is the header file for the micro simulation class.
// It is compiled with nanobind to create a python module, which is then imported by the Micro Manager.

#pragma once
#include <iostream>
#include <vector>
#include <variant>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>     // numpy arrays
#include <nanobind/stl/string.h>  // std::string conversion

#include "general.h"
#include "matmodel.h"
#include "MaterialManager.h"
#include "solver.h"

namespace nb = nanobind;

class MicroSimulation {
  public:
    MicroSimulation(int sim_id, bool late_init = false, const std::string &input_file = "input.json", const std::string &config_file = "pyfans-config.json");
    nb::dict solve(const nb::dict &macro_write_data, double dt);

    nb::dict get_state();
    void     set_state(const nb::dict &state);

    int  get_id();
    void set_id(int id);

  private:
    int         _sim_id;
    std::string _input_file;
    Reader      reader;
    // Hardcoding mechanical models because these definitions need information from the input file.
    using matmanager_t = std::variant<MaterialManager<3, 6> *, MaterialManager<3, 9> *>;
    using solver_t     = std::variant<Solver<3, 6> *, Solver<3, 9> *>;
    matmanager_t matmanager;
    solver_t     solver;
    double       pert_param = 1e-6; // scalar strain perturbation parameter
};

struct PyFANSConfig {
    bool        disable_mpi = false;
    Log::Config logging;
};
