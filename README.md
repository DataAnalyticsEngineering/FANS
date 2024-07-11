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

### Docker
We provide a set of docker images for different use cases:
- fans-ci: Contains the minimum tools to build FANS (including dev packages of dependencies with the required headers), but does not include FANS itself. Meant for a CI workflow.
- fans-dev: Based upon fans-ci, but offers a non-root user and handling of UID and GID to not mess up permissions when volume mounting into the container. Meant for developers that can't install the required tools on their machines.

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
    cmake .
    ```

3. Compile the project using `make`:
    ```sh
    make all
    ```

The compilation will generate a binary called `FANS` in the `test/` directory.

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



