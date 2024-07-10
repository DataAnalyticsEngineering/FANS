set(CMAKE_MODULE_PATH_save "${CMAKE_MODULE_PATH}")
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_LIST_DIR}/customFindScripts")

if ("$ENV{SETVARS_COMPLETED}" STREQUAL "1")
    message(
            WARNING
            "Intel OneAPI environment is active. Might lead to issues with correct MPI lib discovery!"
    )
endif ()

include(CMakeFindDependencyMacro)
find_dependency(HDF5 COMPONENTS CXX)
find_dependency(Eigen3)
find_dependency(OpenMP)
find_dependency(MPI)
find_dependency(FFTW3 COMPONENTS SINGLE DOUBLE LONGDOUBLE MPI)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH_save}")
unset(CMAKE_MODULE_PATH_save)

include(${CMAKE_CURRENT_LIST_DIR}/FANS.cmake)
