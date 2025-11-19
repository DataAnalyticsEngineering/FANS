<p align="center">
  <img src="docs/images/FANS_logo.png" alt="Example Image" width="250" height="250">
</p>

<p align="center">
  <a href="https://prefix.dev/channels/conda-forge/packages/fans"><img src="https://img.shields.io/github/v/release/DataAnalyticsEngineering/FANS?label=Release&color=004191" alt="GitHub Release"></a>
  <a href="https://anaconda.org/conda-forge/fans"><img src="https://anaconda.org/conda-forge/fans/badges/platforms.svg" alt="Anaconda-Server Badge"></a>
  <a href="https://github.com/DataAnalyticsEngineering/FANS/actions"><img src="https://github.com/DataAnalyticsEngineering/FANS/workflows/Build%20and%20test%20pixi-build/badge.svg" alt="Build and test pixi-build"></a>
  <a href="https://anaconda.org/conda-forge/fans"><img src="https://anaconda.org/conda-forge/fans/badges/downloads.svg" alt="Anaconda-Server Badge"></a>
  <a href="https://prefix.dev/channels/conda-forge/packages/fans"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json" alt="Pixi Badge"></a>
  <img src="https://img.shields.io/github/last-commit/DataAnalyticsEngineering/FANS" alt="Last commit">
  <a href="https://github.com/DataAnalyticsEngineering/FANS/blob/main/LICENSE"><img src="https://img.shields.io/github/license/DataAnalyticsEngineering/FANS" alt="License"></a>
  <a href="https://github.com/DataAnalyticsEngineering/FANS/stargazers"><img src="https://img.shields.io/github/stars/DataAnalyticsEngineering/FANS?style=social" alt="Stars"></a>
</p>

# Fourier-Accelerated Nodal Solver (FANS)

Fourier-Accelerated Nodal Solver (FANS) is an FFT-based homogenization solver for microscale multiphysics problems. FANS is written in C++, built using CMake, and it has MPI parallelization.

## Table of Contents

- [Quick start](#quick-start)
- [Build from source](#build-from-source)
  - [Installing dependencies](#installing-dependencies)
  - [Building FANS](#building-fans)
- [Python environment for the FANS dashboard](#python-environment-for-the-fans-dashboard)
- [Input file format](#input-file-format)
  - [Microstructure definition](#microstructure-definition)
  - [Problem type and material model](#problem-type-and-material-model)
  - [Solver settings](#solver-settings)
  - [Macroscale loading conditions](#macroscale-loading-conditions)
  - [Results specification](#results-specification)

## Quick start

**Want to get started immediately?**

FANS is available as a precompiled binary on [conda-forge](https://anaconda.org/conda-forge/fans). Package managers such as [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html), [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html), [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html), and [Pixi](https://pixi.sh) can be used to install FANS from the conda-forge channel.

Use [Pixi](https://pixi.sh) (recommended):

```bash
# Install Pixi (if not already installed)
curl -fsSL https://pixi.sh/install.sh | sh

# Install FANS via Pixi
pixi global install fans

# Verify installation
FANS --version
```

That's it! No dependencies to install, no compilation needed 🚀

To get started immediately, we include ready to use example [input files](test/input_files/) and [microstructures](test/microstructures/) you can use as templates to create your own.

---

## Build from source

**Recommended for:** Developers, contributors, HPC users, or those needing custom builds.

FANS requires the following dependencies:

| Dependency | Purpose | |
| ------------ | --------- | ------------------ |
| **C++ Compiler** | (GCC, Clang, etc.) | C++17 or newer |
| **CMake** | Build system | ≥ 3.21 |
| **MPI** | Parallel computing | (OpenMPI, MPICH, Intel MPI) |
| **HDF5** | Data I/O | **with MPI support** |
| **FFTW3** | FFT computations | **with MPI support** |
| **Eigen3** | Linear algebra | ≥ 3.4 |
| **nlohmann-json** | JSON parsing | ≥ 3.11 |

### Installing dependencies

<details>

<summary><b> Using Pixi (Cross-platform - Easiest for source builds)</b></summary>

This uses the repository's [pixi.toml](pixi.toml) to define the `dev` environment.

```bash
# Clone the repository
git clone https://github.com/DataAnalyticsEngineering/FANS.git
cd FANS

# Enter development environment (all dependencies pre-installed!)
pixi shell -e dev
```

</details>

<details>
<summary><b>Linux (Debian/Ubuntu)</b></summary>

We recommend installing the dependencies using `apt`:

```bash
apt-get install -y \
    build-essential \
    cmake \
    git \
    file \
    libhdf5-dev \
    libhdf5-openmpi-dev \
    libopenmpi-dev \
    libeigen3-dev \
    libfftw3-dev \
    libfftw3-mpi-dev \
    nlohmann-json3-dev
```

</details>

<details>
<summary><b>macOS</b></summary>

We recommend installing the dependencies using [`brew`](https://brew.sh):

```bash
brew install gnu-time cmake gcc@15
brew install open-mpi --build-from-source --cc=gcc-15
brew install hdf5-mpi --build-from-source --cc=gcc-15
brew install fftw eigen nlohmann-json

# Set environment variables
export CC=gcc-15 CXX=g++-15 MPICC=mpicc MPICXX=mpicxx
```

</details>

<details>
<summary><b> Using Spack (HPC environments)</b></summary>

[Spack](https://spack.readthedocs.io/en/latest/) is a flexible package manager for building and managing software stacks in high-performance computing environments. Install Spack by following these [installation instructions](https://spack.readthedocs.io/en/latest/getting_started.html). Once Spack is set up, install the required dependencies:

```bash
spack install cmake
spack install mpi
spack install hdf5+cxx+mpi
spack install eigen
spack install fftw+mpi
spack install nlohmann-json

# Load dependencies
spack load cmake mpi hdf5 eigen fftw nlohmann-json
```

Additionally, optimized FFTW implementations can be used depending on your system's architecture:

- AMD systems: `spack install amdfftw+mpi`
- Cray systems: `spack install cray-fftw+mpi`
- Fujitsu systems: `spack install fujitsu-fftw+mpi`

</details>

<details>
<summary><b>Docker images</b></summary>

Pre-configured Docker images are available for containerized deployments. See [`docker/README.md`](docker/README.md) for further details.

</details>

### Building FANS

```bash
# Clone the repository
git clone https://github.com/DataAnalyticsEngineering/FANS.git
cd FANS

# Create build directory
mkdir build && cd build

# Configure (basic)
cmake ..

# Build
cmake --build . -j

# Run tests with 8 mpi processes
cd ../test
./run_tests.sh -n 8
```

**Build options:**

| CMake Option | Description | Default |
| -------------- | ------------- | --------- |
| `CMAKE_BUILD_TYPE` | Build type: `Debug`, `Release`, `RelWithDebInfo` | `NONE` |
| `CMAKE_INTERPROCEDURAL_OPTIMIZATION` | Enable link-time optimization (LTO) | `ON` (if supported) |
| `FANS_BUILD_STATIC` | Build static library | `OFF` |
| `CMAKE_INSTALL_PREFIX` | Installation directory | System default |
| `FANS_LIBRARY_FOR_MICRO_MANAGER` | Build Python bindings using Pybind11 (needed) | `OFF` |
| `FANS_ENABLE_SANITIZERS` | Enable runtime sanitizers (AddressSanitizer and LeakSanitizer) for memory debugging | `OFF` |

---

## Python environment for the FANS dashboard

FANS includes [`FANS_Dashboard.ipynb`](FANS_Dashboard/FANS_Dashboard.ipynb), a comprehensive pipeline for post-processing, visualization, and analysis of simulation results. We recommend setting up a Python virtual environment via [Pixi](https://pixi.sh/) with all required Python dependencies in an isolated environment:

```bash
# Install and activate the dashboard environment
pixi shell -e dashboard
```

The `dashboard` environment includes:

- Python
- Jupyter notebook (`ipykernel`)
- [MSUtils](https://github.com/DataAnalyticsEngineering/MSUtils) for FANS-specific utilities
- Testing tools (`pytest`)
- Code quality tools (`pre-commit`)

See [`FANS_Dashboard`](FANS_Dashboard/) for further details.

---

## Input file format

FANS requires a JSON input file specifying the problem parameters. Example input files can be found in the [`test/input_files`](test/input_files) directory. It is recommended to use these files as a reference to create your input file.

### Microstructure definition

```json
"microstructure": {
                    "filepath": "microstructures/sphere32.h5",
                    "datasetname": "/sphere/32x32x32/ms",
                    "L": [1.0, 1.0, 1.0]
                  }
```

- `filepath`: This specifies the path to the HDF5 file that contains the microstructure data.
- `datasetname`: This is the path within the HDF5 file to the specific dataset that represents the microstructure.
- `L`: Microstructure length defines the physical dimensions of the microstructure in the $x$, $y$, and $z$ directions.

### Problem type and material model

```json
"problem_type": "mechanical",
"strain_type": "small",
"materials": [
    {
        "phases": [0],
        "matmodel": "PseudoPlasticLinearHardening",
        "material_properties": {
            "bulk_modulus": [62.5000],
            "shear_modulus": [28.8462],
            "yield_stress": [0.1],
            "hardening_parameter": [0.0]
        }
    },
    {
        "phases": [1],
        "matmodel": "LinearElasticIsotropic",
        "material_properties": {
            "bulk_modulus": [222.222],
            "shear_modulus": [166.6667]
        }
    }
]
```

- `problem_type`: This defines the type of physical problem you are solving. Options include `thermal` problems and `mechanical` problems.
- `strain_type`: This indicates whether the problem is formulated using infinitesimal (`small`) strain or finite (`large`) strain theory.
- `materials`: An array of material groups, where each group assigns one or more phases to a specific material model. Each material group contains:
  - `phases`: An array of phase IDs (material labels) from the microstructure that use this material model.
  - `matmodel`: The constitutive model for this material group. Available models include:

    - `LinearThermalIsotropic` for linear isotropic conductive material model.
    - `LinearThermalTriclinic` for linear triclinic conductive material model.
    - `GBDiffusion` for diffusion model with transversely isotropic grain boundary and isotropic bulk for polycrystalline materials.

    - `LinearElasticIsotropic` for linear isotropic elastic material model.
    - `LinearElasticTriclinic` for linear triclinic elastic material model.
    - `PseudoPlasticLinearHardening` / `PseudoPlasticNonLinearHardening` for plasticity mimicking model with linear/nonlinear hardening.
    - `J2ViscoPlastic_LinearIsotropicHardening` / `J2ViscoPlastic_NonLinearIsotropicHardening` for rate-independent / dependent J2 plasticity model with kinematic and linear/nonlinear isotropic hardening.
    - `SaintVenantKirchhoff` for the hyperelastic Saint Venant-Kirchhoff material model.
    - `CompressibleNeoHookean` for the compressible Neo-Hookean material model.

  - `material_properties`: Material parameters specific to the chosen model. Properties are defined as arrays, where each element corresponds to one of the phases listed in the `phases` array.

### Solver settings

```json
"FE_type": "HEX8",
"method": "cg",
"error_parameters":{
                     "measure": "Linfinity",
                     "type": "absolute",
                     "tolerance": 1e-10
                   },
"n_it": 100,
```

- `FE_type`: This specifies the type of finite element to be used. Common options include:
  - `HEX8`: Standard trilinear hexahedral elements with full integration (8 Gauss points). Suitable for most problems but may exhibit volumetric locking for nearly incompressible materials (Poisson's ratio ~ 0.5).
  - `BBAR`: B-bar elements with selective reduced integration to mitigate volumetric locking. Recommended for materials with high Poisson's ratios (0.4 to 0.5).
  - `HEX8R`: Reduced integration elements with a single Gauss point at the element center. Use with caution—these may produce less accurate field results and can cause local material issues such as negative Jacobian ($J < 0$), leading to nonphysical solutions (hourglassing).
- `method`: This indicates the numerical method to be used for solving the system of equations. `cg` stands for the Conjugate Gradient method, and `fp` stands for the Fixed Point method.
- `error_parameters`: This section defines the error parameters for the solver. Error control is applied to the finite element nodal residual of the problem.
  - `measure`: Specifies the norm used to measure the error. Options include `Linfinity`, `L1`, or `L2`.
  - `type`: Defines the type of error measurement. Options are `absolute` or `relative`.
  - `tolerance`: Sets the tolerance level for the solver, defining the convergence criterion based on the chosen error measure. The solver iterates until the solution meets this tolerance.
- `n_it`: Specifies the maximum number of iterations allowed for the FANS solver.

### Macroscale loading conditions

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

- `macroscale_loading`: This defines the external loading applied to the microstructure. It is an array of arrays, where each sub-array represents a load path applied to the system. The format of the load path depends on the problem type:
  - For `thermal` problems, the array typically has 3 components, representing the macroscale temperature gradients in the $x$, $y$, and $z$ directions.
  - For `small` strain `mechanical` problems, the array must have 6 components, corresponding to the macroscale strain tensor in Mandel notation: $[\varepsilon_{11},\ \varepsilon_{22},\ \varepsilon_{33},\ \sqrt{2}\,\varepsilon_{12},\ \sqrt{2}\,\varepsilon_{13},\ \sqrt{2}\,\varepsilon_{23}]$.
  - For `large` strain `mechanical` problems, the array must have 9 components, corresponding to the macroscale deformation gradient tensor: $[F_{11},\ F_{12},\ F_{13},\ F_{21},\ F_{22},\ F_{23},\ F_{31},\ F_{32},\ F_{33}]$.

In the case of path/time-dependent loading, as shown, for example, in plasticity problems, the `macroscale_loading` array can include multiple steps with corresponding loading conditions.

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

### Results specification

```json
"results": ["stress_average", "strain_average", "absolute_error", "phase_stress_average", "phase_strain_average",
            "microstructure", "displacement", "displacement_fluctuation", "stress", "strain"]
```

- `results`: This array lists the quantities that should be stored in the results HDF5 file during the simulation. Each string in the array corresponds to a specific result:

  - `stress_average` and `strain_average`: Volume averaged- homogenized stress and strain over the entire microstructure.
  - `absolute_error`: The L-infinity error of the finite element nodal residual at each iteration.
  - `phase_stress_average` and `phase_strain_average`: Volume averaged- homogenized stress and strain for each phase within the microstructure.
  - `microstructure`: The original microstructure data.
  - `displacement`: The displacement field (for mechanical problems) and temperature field (for thermal problems) at each voxel in the microstructure.
  - `displacement_fluctuation`: The periodic displacement fluctuation field (for mechanical problems) and periodic temperature fluctuation field (for thermal problems at each voxel in the microstructure).
  - `stress` and `strain`: The stress and strain fields at each voxel in the microstructure.

- Additional material model-specific results can be included depending on the problem type and material model.

## Acknowledgements

Funded by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany’s Excellence Strategy - EXC 2075 – 390740016. Contributions by Felix Fritzen are funded by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) within the Heisenberg program - DFG-FR2702/8 - 406068690; DFG-FR2702/10 - 517847245 and through NFDI-MatWerk - NFDI 38/1 - 460247524. We acknowledge the support of the Stuttgart Center for Simulation Science ([SimTech](https://www.simtech.uni-stuttgart.de/)).

## Contributors

- [Sanath Keshav](https://github.com/sanathkeshav)
- [Florian Rieg](https://github.com/scylent)
- [Ishaan Desai](https://github.com/IshaanDesai)
- [Moritz Sigg](https://github.com/siggmo)
- [Claudius Haag](https://github.com/claudiushaag)
- [Felix Fritzen](https://github.com/EMMAOpenSource)
