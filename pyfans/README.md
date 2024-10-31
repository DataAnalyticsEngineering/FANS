# pyFANS

pyFANS is a Python-wrapped library to control FANS via the [Micro Manager](https://precice.org/tooling-micro-manager-overview.html). The main idea is to create a large number of FANS simulations, and couple them to one macro-scale simulation typically in Abaqus, CalculiX, etc. The library follows the [API of the Micro Manager](https://precice.org/tooling-micro-manager-prepare-micro-simulation.html).

## Dependencies

- [pybind11](https://pybind11.readthedocs.io/en/stable/index.html)

## Building

To build FANS as a library, set the CMake variable `FANS_LIB` to `ON`. The CMake command to compile FANS would then be `cmake .. -DFANS_LIB=ON`.

## Usage

pyFANS is intended to be used with the Micro Manager and preCICE for two-scale coupled simulations. However, standalone use of the library is not restricted per se. Look at the [test_pyfans](../test/test_pyfans/) example to see how the library is used in a Python script.
