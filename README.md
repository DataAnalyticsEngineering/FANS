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

### Docker
We provide a set of docker images to ease the installation process:
1. `unistuttgartdae/fans-runtime`
   - **Description**: This image provides the runtime environment for FANS, including all necessary libraries and dependencies required to execute FANS. It is based on Ubuntu 22.04 and includes libraries like HDF5, OpenMPI, Eigen3, and FFTW3.

2. `unistuttgartdae/fans-build-env`
   - **Description**: This image includes the build environment for FANS. It contains all tools and libraries needed to compile the FANS source code, including build essentials, CMake, and Git. It builds upon the `fans-runtime` image.

3. `unistuttgartdae/fans-build`
   - **Description**: This image contains the built FANS application. It copies the entire FANS repository into the image, compiles the code using CMake, and packages it. This image is used internally to create the final executable package.

4. `unistuttgartdae/fans-publish`
   - **Description**: This image is the final runtime environment for FANS with the built package installed. It is derived from the `fans-runtime` image, installs the compiled FANS package, and sets up a non-root user environment. This image is intended for end-users to run the FANS application.

5. `unistuttgartdae/fans-dev-env`
   - **Description**: This image provides a development environment for FANS. It also sets up a non-root user environment for developers. This image is based on the `fans-build-env` image and is intended for development and debugging purposes.

Only the last two images are intended to actually be used. To set up a container based upon them, type:
```bash
mkdir workspace
cd workspace
[git clone -b CMake git@github.tik.uni-stuttgart.de:DAE/FANS.git]

docker pull unistuttgartdae/fans-[publish/dev-env]
docker create --name fans -it \
  -e HOST_UID=$(id -u) \
  -e HOST_GID=$(id -g) \
  -v /etc/localtime:/etc/localtime:ro \
  -v /etc/timezone:/etc/timezone:ro \
  -v $PWD/:/workspace/ \
  unistuttgartdae/fans-[publish/dev-env]
```

Then, to start the container and attach a shell, run:
```bash
docker start -i fans
```

In case of the image `fans-publish`, to directly execute FANS without attaching to the container, use (be aware that the paths must be relative to the container!):
```bash
docker run fans [arguments]
```

In case of the image `fans-dev-env`, the intended workflow would be to edit the code as usual in your IDE, but then to build it you would attach to the container and run:
```bash
cd FANS
mkdir build
cd build
cmake ..
cmake --build . -j

```

## Usage

To run the FANS solver, you need to provide a JSON input file specifying the problem parameters. Here is the command to run FANS:

```sh
nohup /usr/bin/time -v mpiexec -n 1 ./FANS path/to/your/input_file.json &
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



