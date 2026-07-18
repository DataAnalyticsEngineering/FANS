# pyFANS

pyFANS is a Python-wrapped library to control FANS via the [Micro Manager](https://precice.org/tooling-micro-manager-overview.html). The main idea is to create a large number of FANS simulations, and couple them to one macro-scale simulation typically in Abaqus, CalculiX, etc. The library follows the [API of the Micro Manager](https://precice.org/tooling-micro-manager-prepare-micro-simulation.html).

## Dependencies

- [pybind11](https://pybind11.readthedocs.io/en/stable/index.html)

## Building

To build FANS as a Micro Manager compatible Python library, set the CMake variable `FANS_LIB` to `ON`. The CMake command to compile FANS would then be `cmake .. -DFANS_LIBRARY_FOR_MICRO_MANAGER=ON`.

To build multiple pyFANS targets configured to load different input and configuration files, respective build directories must be created. Within those, the corresponding CMAKE flags must be set.
In the following example we show how to configure two different builds:

__NOTE__: `MicroSimulationSuffix` may only be set to non-negative integer values within 0 to 10

In the project root directory, create two build folders:

```bash
cd <FANS-ROOT>
mkdir build0
mkdir build1
```

To configure and build the first target:

```bash
cd build0
cmake -DFANS_LIBRARY_FOR_MICRO_MANAGER=ON -DMicroSimulationSuffix=0 ..
cmake --build -j <N> ..
cp pyfans/PyFANS.*.so <Simulation Directory>/PyFANS0.so
cd ..
```

This will result in the following configuration:

| Name                 | Value               |
|----------------------|---------------------|
| CLASS_NAME           | MicroSimulation0    |
| INPUT_FILE           | input0.json         |
| CONFIG_FILE          | pyfans-config0.json |
| pybind11 MODULE_NAME | PyFANS0             |

To configure and build the second target:

```bash
cd build1
cmake -DFANS_LIBRARY_FOR_MICRO_MANAGER=ON -DMicroSimulationSuffix=1 ..
cmake --build -j <N> ..
cp pyfans/PyFANS.*.so <Simulation Directory>/PyFANS1.so
cd ..
```

This will result in the following configuration:

| Name                 | Value               |
|----------------------|---------------------|
| CLASS_NAME           | MicroSimulation1    |
| INPUT_FILE           | input1.json         |
| CONFIG_FILE          | pyfans-config1.json |
| pybind11 MODULE_NAME | PyFANS1             |

## Usage

pyFANS is intended to be used with the Micro Manager and preCICE for two-scale coupled simulations. However, standalone use of the library is not restricted per se. Look at the [test_pyfans](../test/test_pyfans/) example to see how the library is used in a Python script.

To disable parallelization with MPI, provide a `pyfans-config.json` file containing the field: `no_mpi: true`.
It should be located in the same directed from which for example the Micro Manager is launched.
