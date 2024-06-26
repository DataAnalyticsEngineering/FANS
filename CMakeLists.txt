cmake_minimum_required(VERSION 3.0)
project(FANS)

set(CMAKE_C_COMPILER   mpicc)
set(CMAKE_CXX_COMPILER mpic++)

# with -fno-implicit-templates I get linker errors when using std:: stuff
set(CMAKE_CXX_FLAGS    "-O3 -mavx2 -mfma -fopenmp -lm -DH5_HAVE_PARALLEL")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY test/)
set(HDF5_ENABLE_PARALLEL ON)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package (Eigen3 REQUIRED)


include_directories(${HDF5_INCLUDE_DIRS}
                    /usr/local/include
                    /usr/include
                    include/)

link_directories(/usr/lib/x86_64-linux-gnu/)

add_executable (FANS
                src/main.cpp
                src/matmodel.cpp
                src/reader.cpp)

target_link_libraries(FANS ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} fftw3 fftw3_mpi Eigen3::Eigen)
