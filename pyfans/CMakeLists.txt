# pyFANS builds against an installed FANS, so the FANS core needs no Python.

cmake_minimum_required(VERSION 3.21)

# C too: FANS's FindFFTW3 probes with a try_run on a .c file.
project(pyfans LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(FANS REQUIRED)

# Development.SABIModule provides Python::SABIModule; without it nanobind
# silently ignores STABLE_ABI and emits a version-specific module.
find_package(Python 3.12 REQUIRED
             COMPONENTS Interpreter Development.Module Development.SABIModule)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(PyFANS STABLE_ABI NB_STATIC micro.hpp micro.cpp)
target_link_libraries(PyFANS PRIVATE FANS::FANS)
target_include_directories(PyFANS PRIVATE ${HDF5_INCLUDE_DIRS} ${FFTW3_INCLUDE_DIRS})

# From site-packages, libFANS.so is two levels up in <prefix>/lib.
if (APPLE)
    set_target_properties(PyFANS PROPERTIES INSTALL_RPATH "@loader_path/../..")
else ()
    set_target_properties(PyFANS PROPERTIES INSTALL_RPATH "$ORIGIN/../..")
endif ()

install(TARGETS PyFANS LIBRARY DESTINATION .)
