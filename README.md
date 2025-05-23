# Fourier-Accelerated Nodal Solvers (FANS)

Fourier-Accelerated Nodal Solver (FANS) is an FFT-based homogenization solver for microscale multiphysics problems. FANS is written in C++, built using CMake, and it has MPI parallelization.

<p align="center">
  <img src="docs/images/FANS_logo.png" alt="Example Image" width="400" height="400">
</p>

## Table of contents

- [Dependencies](#dependencies)
- [Building](#building)
- [Installing](#installing)
- [Input File Format](#input-file-format)

## Dependencies

FANS has the following dependencies:

- A C++ compiler (e.g. GCC, Clang, etc.)
- CMake (version 3.21 or higher)
- Git (for cloning this repo)
- MPI (mpicc and mpic++)
- HDF5 with MPI support
- Eigen3
- FFTW3 with MPI support

### Installing dependencies

- On Debian based systems, we recommend installing the dependencies using using `apt`,

  ```bash
  apt-get install \
      libhdf5-dev \
      libopenmpi-dev \
      libeigen3-dev \
      libfftw3-dev \
      libfftw3-mpi-dev
  ```

- On macOS, you can obtain the dependencies using `brew` and set the environment variables:

  ```zsh
  brew install gnu-time cmake gcc@14
  brew install open-mpi --build-from-source --cc=gcc-14
  brew install hdf5-mpi --build-from-source --cc=gcc-14
  brew install fftw eigen

  export CC=gcc-14 CXX=g++-14 MPICC=mpicc MPICXX=mpicxx
  ```

### Setting up a Python environment

Also, we recommend to set up a Python virtual environment for the [`FANS_Dashboard.ipynb`](FANS_Dashboard/FANS_Dashboard.ipynb) via [pixi](https://pixi.sh/) with all required Python dependencies in an isolated environment:

```bash
# Install pixi if not done already
curl -fsSL https://pixi.sh/install.sh | sh

# Create and activate the environment
pixi shell
```

We also provide a set of Docker images. For further information, please refer to the [Docker README](docker/README.md).

### Installing dependencies using Spack

Spack is a package manager designed for high-performance computing environments. It simplifies the installation of complex software stacks, making it ideal for setting up FANS on HPC systems.

1. **Install Spack** by following these [installation instructions](https://spack.readthedocs.io/en/latest/getting_started.html).

2. **Install Dependencies**: Once Spack is set up, install the required dependencies:

    ```bash
    spack install cmake
    spack install mpi
    spack install hdf5 +cxx +mpi
    spack install eigen
    spack install fftw +mpi
    ```

    Additionally, optimized FFTW implementations can be used depending on your system's architecture, for example `amdfftw` (For AMD systems) or `cray-fftw` (For Cray systems) or `fujitsu-fftw` (For Fujitsu systems).

3. **Load Dependencies** Once dependencies are installed, load them before building:

    ```bash
    spack load cmake mpi hdf5 eigen fftw
    ```

## Building

1. Clone the repository:

    ```bash
    git clone https://github.com/DataAnalyticsEngineering/FANS.git
    cd FANS
    ```

2. Configure the build using CMake:

    ```bash
    mkdir build
    cd build
    cmake ..
    ```

3. Compile:

    ```bash
    cmake --build . -j
    ```

The compilation symlinks the generated `FANS` binary into the `test/` directory for convenience.

### Configuring a build

The following CMake configuration options exist:

- `CMAKE_BUILD_TYPE`: Sets the build type. Common values are Debug, Release, RelWithDebInfo, and MinSizeRel.
  - Default: NONE

- `FANS_BUILD_STATIC`: Build static library instead of shared library.
  - Default: OFF

- `CMAKE_INTERPROCEDURAL_OPTIMIZATION`: Enable inter-procedural optimization (IPO) for all targets.
  - Default: ON (if supported)
  - Note: When you run the configure step for the first time, IPO support is automatically checked and enabled if available. A status message will indicate whether IPO is activated or not supported.

## Installing

Install FANS (system-wide) using the following options:

1. Using CMake (sudo required if --prefix is omitted):

    ```bash
    cmake --install . [--prefix <install-dir>]
    ```

### Install using Conda

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/fans/badges/version.svg)](https://anaconda.org/conda-forge/fans)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/fans/badges/platforms.svg)](https://anaconda.org/conda-forge/fans)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/fans/badges/downloads.svg)](https://anaconda.org/conda-forge/fans)

FANS is also available as a conda-package on [conda-forge/fans](https://anaconda.org/conda-forge/fans). No dependencies have to be manually installed for it to work.
It can be installed via

```bash
conda install conda-forge::fans
```

exposing the executable `FANS`.

## Input File Format

FANS requires a JSON input file specifying the problem parameters. Example input files can be found in the [`test/input_files`](test/input_files) directory. It is recommended to use these files as a reference to create your own input file.

### Microstructure Definition

```json
"microstructure": {
                    "filepath": "microstructures/sphere32.h5",
                    "datasetname": "/sphere/32x32x32/ms",
                    "L": [1.0, 1.0, 1.0]
                  }
```

- `filepath`: This specifies the path to the HDF5 file that contains the microstructure data.
- `datasetname`: This is the path within the HDF5 file to the specific dataset that represents the microstructure.
- `L`: Microstructure length defines the physical dimensions of the microstructure in the x, y, and z directions.

### Problem Type and Material Model

```json
"matmodel": "LinearElasticIsotropic",
"material_properties": {
                         "bulk_modulus": [62.5000, 222.222],
                         "shear_modulus": [28.8462, 166.6667]
                       }
```

- `problem_type`: This defines the type of physical problem you are solving. Common options include `thermal` problems and `mechanical` problems.
- `matmodel`: This specifies the material model to be used in the simulation. Examples include

  - `LinearThermalIsotropic` for linear isotropic conductive material model
  - `LinearThermalTriclinic` for linear triclinic conductive material model
  - `GBDiffusion` for diffusion model with transversely isotropic grain boundary and isotropic bulk for polycrystalline materials

  - `LinearElasticIsotropic` for linear isotropic elastic material model
  - `LinearElasticTriclinic` for linear triclinic elastic material model
  - `PseudoPlasticLinearHardening` / `PseudoPlasticNonLinearHardening` for plasticity mimicking model with linear/nonlinear hardening
  - `J2ViscoPlastic_LinearIsotropicHardening` / `J2ViscoPlastic_NonLinearIsotropicHardening` for rate independent / dependent J2 plasticity model with kinematic and linear/nonlinear isotropic hardening.

- `material_properties`: This provides the necessary material parameters for the chosen material model. For thermal problems, you might specify `conductivity`, while mechanical problems might require `bulk_modulus`, `shear_modulus`, and more properties for advanced material models. These properties can be defined as arrays to represent multiple phases within the microstructure.

### Solver Settings

```json
"method": "cg",
"error_parameters":{
                     "measure": "Linfinity",
                     "type": "absolute",
                     "tolerance": 1e-10
                   },
"n_it": 100,
```

- `method`: This indicates the numerical method to be used for solving the system of equations. `cg` stands for the Conjugate Gradient method, and `fp` stands for the Fixed Point method.
- `error_parameters`: This section defines the error parameters for the solver. Error control is applied on the finite element nodal residual of the problem.
  - `measure`: Specifies the norm used to measure the error. Options include `Linfinity`, `L1`, or `L2`.
  - `type`: Defines the type of error measurement. Options are `absolute` or `relative`.
  - `tolerance`: Sets the tolerance level for the solver, defining the convergence criterion based on the chosen error measure. The solver iterates until the solution meets this tolerance.
- `n_it`: Specifies the maximum number of iterations allowed for the FANS solver.

### Macroscale Loading Conditions

```json
"macroscale_loading":   [
                            [
                                [0.004, -0.002, -0.002, 0, 0, 0],
                                [0.008, -0.004, -0.004, 0, 0, 0],
                                [0.012, -0.006, -0.006, 0, 0, 0],
                                [0.016, -0.008, -0.008, 0, 0, 0],
                            ],
                            [
                                [0, 0, 0, 0.002, 0, 0],
                                [0, 0, 0, 0.004, 0, 0],
                                [0, 0, 0, 0.006, 0, 0],
                                [0, 0, 0, 0.008, 0, 0],
                            ]
                        ],
```

- `macroscale_loading`: This defines the external loading applied to the microstructure. It is an array of arrays, where each sub-array represents a loading condition applied to the system. The format of the loading array depends on the problem type:
- For `thermal` problems, the array typically has 3 components, representing the temperature gradients in the $x$, $y$, and $z$ directions.
- For `mechanical` problems, the array must have 6 components, corresponding to the components of the strain tensor in Mandel notation (e.g., $[\varepsilon_{11}, \varepsilon_{22}, \varepsilon_{33}, \sqrt{2}\varepsilon_{12}, \sqrt{2}\varepsilon_{13}, \sqrt{2}\varepsilon_{23}]$).

In the case of path/time-dependent loading as shown, for example as in plasticity problems, the `macroscale_loading` array can include multiple steps with corresponding loading conditions.

FANS also supports mixed boundary conditions, where some components can be strain-controlled while others are stress-controlled:

```json
"macroscale_loading":   [{
                           "strain_indices" : [2,3,4,5],
                           "stress_indices" : [0,1],
                           "strain" : [[0.005 , 0.0, 0.0, 0.0],
                                       [0.010 , 0.0, 0.0, 0.0]],
                           "stress" : [[0.0, 0.0],
                                       [0.0, 0.0]]
                          }]
```

### Results Specification

```json
"results": ["stress_average", "strain_average", "absolute_error", "phase_stress_average", "phase_strain_average",
            "microstructure", "displacement", "displacement_fluctuation", "stress", "strain"]
```

- `results`: This array lists the quantities that should be stored into the results HDF5 file during the simulation. Each string in the array corresponds to a specific result:

  - `stress_average` and `strain_average`: Volume averaged- homogenized stress and strain over the entire microstructure.
  - `absolute_error`: The L-infinity error of finite element nodal residual at each iteration.
  - `phase_stress_average` and `phase_strain_average`: Volume averaged- homogenized stress and strain for each phase within the microstructure.
  - `microstructure`: The original microstructure data.
  - `displacement`: The displacement field (for mechanical problems) and temperature field (for thermal problems) at each voxel in the microstructure.
  - `displacement_fluctuation`: The periodic displacement fluctuation field (for mechanical problems) and periodic temperature fluctuation field (for thermal problems at each voxel in the microstructure).
  - `stress` and `strain`: The stress and strain fields at each voxel in the microstructure.

- Additional material model specific results can be included depending on the problem type and material model.

## Acknowledgements

Funded by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany’s Excellence Strategy - EXC 2075 – 390740016. Contributions by Felix Fritzen are funded by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) within the Heisenberg program - DFG-FR2702/8 - 406068690; DFG-FR2702/10 - 517847245 and through NFDI-MatWerk - NFDI 38/1 - 460247524. We acknowledge the support by the Stuttgart Center for Simulation Science ([SimTech](https://www.simtech.uni-stuttgart.de/)).

## Contributors

- [Sanath Keshav](https://github.com/sanathkeshav)
- [Florian Rieg](about:blank)
- [Ishaan Desai](https://github.com/IshaanDesai)
- [Moritz Sigg](https://github.com/siggmo)
- [Claudius Haag](https://github.com/claudiushaag)
- [Felix Fritzen](https://github.com/EMMAOpenSource)
