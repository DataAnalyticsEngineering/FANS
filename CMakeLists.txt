cmake_minimum_required(VERSION 3.0...3.28)

# ##############################################################################
# GENERAL SETTINGS
# ##############################################################################

project(FANS
        VERSION 0.3.0
        LANGUAGES C CXX
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_CXX_STANDARD 14)

# with -fno-implicit-templates I get linker errors when using std:: stuff
# TODO: should be developer specific, by using e.g. CMake Presets
set(CMAKE_CXX_FLAGS "-fmax-errors=1 -O3")
set(CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH})

# IPO
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.9")
    cmake_policy(SET CMP0069 NEW)
    include(CheckIPOSupported)
    get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
    check_ipo_supported(RESULT result LANGUAGES ${languages})
    if (result AND NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION)
        option(CMAKE_INTERPROCEDURAL_OPTIMIZATION
                "Enable interprocedural optimization for all targets." ON
        )
        message(STATUS "IPO activated.")
    elseif (NOT result)
        message(STATUS "IPO not supported.")
    endif ()
endif ()

# ##############################################################################
# RPATH settings
# ##############################################################################

# Taken from https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# From https://stackoverflow.com/questions/73248130/how-to-avoid-the-removal-of-the-rpath-during-cmake-install-step

# Set RPATH to be relative and honor user overrides, whether at the command line or FetchContent
if(APPLE)
        set(RPATH_BASE "@loader_path")
else()
        set(RPATH_BASE "$ORIGIN")
endif()
file(RELATIVE_PATH REL_PATH_LIB
        "/${CMAKE_INSTALL_BINDIR}/"
        "/${CMAKE_INSTALL_LIBDIR}/")
set(CMAKE_INSTALL_RPATH "${RPATH_BASE};${RPATH_BASE}/${REL_PATH_LIB}"
        CACHE STRING "Install RPATH")

# ##############################################################################
# REQUIRED PACKAGES
# ##############################################################################

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

if ("$ENV{SETVARS_COMPLETED}" STREQUAL "1")
    message(
            WARNING
            "Intel OneAPI environment is active, which might lead to issues with MPI discovery."
    )
endif ()

set(HDF5_ENABLE_PARALLEL ON)
find_package(HDF5 REQUIRED COMPONENTS CXX)

find_package(Eigen3 REQUIRED)

find_package(OpenMP REQUIRED)

find_package(MPI REQUIRED)

find_package(FFTW3 REQUIRED COMPONENTS DOUBLE MPI)

option(FANS_LIBRARY_FOR_MICRO_MANAGER "Building FANS as a library to be used by the Micro Manager." OFF)

if (FANS_LIBRARY_FOR_MICRO_MANAGER)
    include(FetchContent)
    FetchContent_Declare(
        pybind11 GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v2.12.0
    )
    FetchContent_MakeAvailable(pybind11)
endif()

# ##############################################################################
# TARGETS
# ##############################################################################

option(FANS_BUILD_STATIC "Build static library" OFF)
if (FANS_BUILD_STATIC)
    add_library(FANS_FANS STATIC)
else ()
    add_library(FANS_FANS SHARED)
endif ()
add_library(FANS::FANS ALIAS FANS_FANS)

if (CMAKE_SYSTEM_PROCESSOR MATCHES "arm|aarch64")
    target_compile_options(FANS_FANS PUBLIC -march=armv8-a+simd+fp+crypto)
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
    # Detecting if there is avx2 and fma support on linux
    # Deteced together, because these flags were introduced by Intel Haskell 4th. Gen
    # Distinguishing between mac/linux because lscpu is used, which is not supported on mac.
    if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        # macOS: Assume AVX2 and FMA are supported as per requirement
        message(STATUS "Assuming AVX2 and FMA support on macOS; using -mavx2 and -mfma.")
        target_compile_options(FANS_FANS PUBLIC -mavx2 -mfma)
    elseif (CMAKE_SYSTEM_NAME STREQUAL "Linux")
        # Linux: Check for AVX2 support using lscpu
        execute_process(COMMAND lscpu OUTPUT_VARIABLE LSCPU_OUTPUT ERROR_QUIET)

        # Default to AVX; enable AVX2 if supported by hardware
        set(FANS_USE_AVX2 FALSE)
        if (LSCPU_OUTPUT MATCHES "avx2")
            set(FANS_USE_AVX2 TRUE)
        endif()

        # Apply compiler options based on hardware detection
        if (FANS_USE_AVX2)
            message(STATUS "AVX2 and FMA supported by hardware; using -mavx2 and -mfma.")
            target_compile_options(FANS_FANS PUBLIC -mavx2 -mfma)
        else()
            message(WARNING "AVX2 and FMA not supported by hardware; falling back to AVX.")
            target_compile_options(FANS_FANS PUBLIC -mavx)
        endif()
    else()
        message(WARNING "Unsupported system for AVX2 detection; defaulting to -mavx.")
        target_compile_options(FANS_FANS PUBLIC -mavx)
    endif()
endif()

add_executable(FANS_main)
add_custom_command(
        TARGET FANS_main
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
        $<TARGET_FILE:FANS_main>
        ${CMAKE_CURRENT_SOURCE_DIR}/test/$<TARGET_FILE_NAME:FANS_main>
        COMMENT "Create a symlink for FANS executable to ${CMAKE_CURRENT_SOURCE_DIR}/test/"
)

if (FANS_LIBRARY_FOR_MICRO_MANAGER)
    add_subdirectory(pyfans)
endif ()

# ##############################################################################
# HEADERS
# ##############################################################################

target_include_directories(FANS_FANS PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# TODO: when creating a shared library in the future, maybe create one single
# FANS.h header that consuming projects can use
set_property(TARGET FANS_FANS PROPERTY PUBLIC_HEADER
        include/general.h
        include/json.hpp
        include/matmodel.h
        include/reader.h
        include/solverCG.h
        include/solverFP.h
        include/solver.h
        include/setup.h

        include/material_models/LinearThermal.h
        include/material_models/LinearElastic.h
        include/material_models/PseudoPlastic.h
        include/material_models/J2Plasticity.h
)

# version.h is generated and added to the build include directory
configure_file(
    "${PROJECT_SOURCE_DIR}/include/version.h.in"
    "${PROJECT_BINARY_DIR}/include/version.h" @ONLY
)
target_include_directories(FANS_main PRIVATE "${PROJECT_BINARY_DIR}/include")

# ##############################################################################
# SOURCES
# ##############################################################################

target_sources(FANS_FANS PRIVATE
        src/reader.cpp
)

target_sources(FANS_main PRIVATE
        src/main.cpp
)

# ##############################################################################
# LINKING
# ##############################################################################

target_link_libraries(FANS_FANS PRIVATE m)

# TODO: when switching to a newer CMake version this can all be done by one
# call to target_link_libraries(FANS_FANS PUBLIC HDF5::HDF5). But CMake 3.16
# does not yet support this.
target_include_directories(FANS_FANS PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(FANS_FANS PUBLIC ${HDF5_CXX_LIBRARIES})
target_compile_definitions(FANS_FANS PUBLIC ${HDF5_DEFINITIONS})
if (HDF5_IS_PARALLEL)
    target_compile_definitions(FANS_FANS PUBLIC H5_HAVE_PARALLEL)
endif ()

target_link_libraries(FANS_FANS PUBLIC MPI::MPI_CXX)

target_include_directories(FANS_FANS PUBLIC ${FFTW3_INCLUDE_DIRS})
target_link_libraries(FANS_FANS PUBLIC ${FFTW3_LIBRARIES})
target_compile_definitions(FANS_FANS PUBLIC ${FFTW3_DEFINITIONS})

target_link_libraries(FANS_FANS PUBLIC Eigen3::Eigen)

target_link_libraries(FANS_FANS PUBLIC OpenMP::OpenMP_CXX)


target_link_libraries(FANS_main PRIVATE FANS::FANS)

# ##############################################################################
# NAMING
# ##############################################################################

set_target_properties(FANS_FANS PROPERTIES OUTPUT_NAME FANS)
set_target_properties(FANS_FANS PROPERTIES EXPORT_NAME FANS)
set_target_properties(FANS_main PROPERTIES OUTPUT_NAME FANS)

# ##############################################################################
# INSTALLING
# ##############################################################################

install(
        TARGETS FANS_FANS
        EXPORT FANS
        LIBRARY
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        COMPONENT FANS_Runtime
        NAMELINK_COMPONENT FANS_Development
        ARCHIVE
        DESTINATION ${CMAKE_INSTALL_LIBDIR}
        COMPONENT FANS_Development
        INCLUDES
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FANS
        PUBLIC_HEADER
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FANS
        COMPONENT FANS_Development
        PRIVATE_HEADER
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FANS
        COMPONENT FANS_Development
)

install(
        TARGETS FANS_main
        RUNTIME
        DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT FANS_Runtime
)

install(
        EXPORT FANS
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FANS
        NAMESPACE FANS::
        FILE FANSTargets.cmake
        COMPONENT FANS_Development
)

install(
        FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/FindFFTW3.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FANS/modules
        COMPONENT FANS_Development
)

write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/FANSConfigVersion.cmake
        COMPATIBILITY SameMajorVersion
)
configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FANSConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/FANSConfig.cmake"
        INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FANS"
        NO_SET_AND_CHECK_MACRO
        NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
install(
        FILES ${CMAKE_CURRENT_BINARY_DIR}/FANSConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/FANSConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FANS
        COMPONENT FANS_Development
)

# ##############################################################################
# TESTING
# ##############################################################################

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(FANS_ENABLE_TESTING_DEFAULT ON)
endif ()
option(FANS_ENABLE_TESTING "Enable testing" ${FANS_ENABLE_TESTING_DEFAULT})
if (FANS_ENABLE_TESTING)
    enable_testing()
    add_subdirectory(test)
endif ()

# ##############################################################################
# PACKAGING
# ##############################################################################

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(FANS_ENABLE_PACKAGING_DEFAULT ON)
endif ()
option(FANS_ENABLE_PACKAGING "Enable creation of Debian packages" ${FANS_ENABLE_PACKAGING_DEFAULT})
if (FANS_ENABLE_PACKAGING)
    add_subdirectory(cmake/packaging)
endif ()
