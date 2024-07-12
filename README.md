# Fourier Accelerated Nodal Solvers (FANS)

Fourier Accelerated Nodal Solvers (FANS) is an FFT-based homogenization solver designed to handle microscale multiphysics problems. This repository contains a C++ implementation of FANS, built using CMake and MPI for parallel computations.

<img src="test/FANS_example.png" alt="Example Image" width="400" height="300">

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Input File Format](#input-file-format)
- [Example](#example)

## Installation

### Prerequisites

Please ensure you have the following dependencies installed on your system:

- CMake (version 3.0 or higher)
- MPI (mpicc and mpic++)
- HDF5 with parallel support
- Eigen3
- FFTW3 with MPI support

Specifically, to run FANS, you at least need the following packages:
```
openmpi-bin, libc6, libfftw3-double3, libfftw3-mpi3, libgcc-s1, libgomp1, libhdf5-103, libopenmpi3, libstdc++6
```
To build fans, you additionally need these packages:
```
libhdf5-dev, libopenmpi-dev, libeigen3-dev, libfftw3-dev, libfftw3-mpi-dev
```

### Docker
We provide a set of docker images for different use cases:
- fans: Contains the minimum environment for FANS to run and has the package 'fans' installed. Offers a non-root user and handling of UID and GID to not mess up permissions when volume mounting into the container. Meant for users of FANS that can't install the fans package directly.
- fans-ci: Contains the minimum tools to build FANS (including dev packages of dependencies with the required headers), but does not include FANS itself. Meant for a CI workflow.
- fans-dev: Based upon fans-ci, but offers a non-root user and handling of UID and GID to not mess up permissions when volume mounting into the container. Meant for developers that can't install the required tools on their machines.

The images are built for both linux/amd64 and linux/arm64 and are available on our (Dockerhub profile)[https://hub.docker.com/u/unistuttgartdae].

To set up a development container with your current working directory mounted into it, type:
```bash
docker create --name fans-dev -it \
  -e HOST_UID=$(id -u) \
  -e HOST_GID=$(id -g) \
  -v /etc/localtime:/etc/localtime:ro \
  -v /etc/timezone:/etc/timezone:ro \
  -v $PWD/:/workspace/ \
  unistuttgartdae/fans-dev
```
The `-e` options provide the entrypoint script of the container with your host user ID and GID, such that the user ID and GID inside the container can be adapted to match yours. This is done to not mess up file permissions in the mounted volumes. The two volume mounts of `/etc/localtime` and `/etc/timezone` are required to have the host date and time inside the container.

To start the container and attach a shell, run:
```bash
docker start -i fans
```
As the `fans-dev` image is meant for developers that can't or don't want to install the required dependencies directly on their machine, the following workflow is suggested: You would work on the code as usual on your host; and only to build and run FANS you would attach to the container (see next section).

### Building the Project

1. Clone the repository:
    ```sh
    git clone https://github.tik.uni-stuttgart.de/SanathKeshav/FANS.git
    cd FANS
    ```
2. Configure the project using CMake:
    ```sh
    mkdir build
    cd build
    cmake ..
    ```

3. Compile the project:
    ```sh
    cmake --build . -j
    ```
The compilation will symlink the generated FANS binary into the `test/` directory for convenience.

### Installing the Project
You can install FANS systemwide using the following options:

1. Using CMake (sudo required if --prefix is omitted):
    ```sh
    cmake --install . [--prefix <install-dir>]
    ```
2. Using .deb packages (only debian based distros; sudo required):
    ```sh
    cpack -G "DEB"
    apt install packages/fans_<version>_<architecture>.deb
    apt install packages/fans-dev_<version>_<architecture>.deb
    ```

### Build Options
This project supports the following CMake build options:

- `CMAKE_BUILD_TYPE`: Sets the build type. Common values are Debug, Release, RelWithDebInfo, and MinSizeRel.

- `FANS_BUILD_STATIC`: Build static library instead of shared library.
  - Default: OFF
  - Usage: `-DFANS_BUILD_STATIC=ON`

- `CMAKE_INTERPROCEDURAL_OPTIMIZATION`: Enable interprocedural optimization (IPO) for all targets.
  - Default: ON (if supported)
  - Usage: `-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF`
  - Note: When you run the configure step for the first time, IPO support is automatically checked and enabled if available. A status message will indicate whether IPO is activated or not supported.

## Usage

To run the FANS solver, you need to provide a JSON input file specifying the problem parameters. Here is the command to run FANS:

```sh
nohup /usr/bin/time -v mpiexec -n 1 ./FANS path/to/your/input_file.json &
```
In case you are using the docker workflow, you first need to start and attach to the container (`docker start -i fans-dev`). If you need to use the command in scripts, the interactive mode (`-i`) is not suitable. Then you need to use `docker exec`:
```bash
docker start fans-dev
nohup /usr/bin/time -v docker exec -u develop -w /workspace/test fans-dev [original command from above] &
docker stop fans-dev
```

## Input File Format

The input file is in JSON format and contains several fields to define the problem settings:
- `comment`: Optional field for comments.
- `ms_filename`: Path to the microstructure file (HDF5 format).
- `ms_datasetname`: Path to the dataset inside the HDF5 file.
- `ms_L`: List defining the domain size.
- `matmodel`: Material model (e.g., `ThermalLinear`, `MechLinear`, `HyperElastic`).
- `material_properties`: Material properties relevant to the chosen material model.
- `problem_type`: Type of the problem (`thermal` or `mechanical`).
- `method`: Solution method (`fp` for fixed point or `cg` for conjugate gradient).
- `TOL`: Tolerance for the solver.
- `n_it`: Maximum number of iterations.
- `g0`: Macroscale loading vector.

## Example

To run a linear elastic mechanical homogenization problem for a single load case on a microstructure image of size `256 x 256 x 256` with a single spherical inclusion,

```sh
mpiexec -n 1 ./FANS input_files/sphere_mech.json
```



