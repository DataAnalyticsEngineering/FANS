cmake_minimum_required(VERSION 3.16)

# ##############################################################################
# GENERAL SETTINGS
# ##############################################################################

project(FANS
        VERSION 0.0.0
        LANGUAGES C CXX
)

# set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# with -fno-implicit-templates I get linker errors when using std:: stuff
# TODO: should be developer specific, by using e.g. CMake Presets
set(CMAKE_CXX_FLAGS "-fmax-errors=1")

execute_process(COMMAND pip show pybind11
                COMMAND grep Location 
                COMMAND cut -d " " -f2
                OUTPUT_VARIABLE SITE_PACKAGES_DIR 
                OUTPUT_STRIP_TRAILING_WHITESPACE)
set(pybind11_DIR "${SITE_PACKAGES_DIR}/pybind11/share/cmake/pybind11/")
message(STATUS "pybind11_DIR: ${pybind11_DIR}")


# Define paths
set(pybind11_INCLUDE_DIR "${SITE_PACKAGES_DIR}/pybind11/include")

# Print the directories for debugging purposes

# Set CMake prefix path to find pybind11Config.cmake
# set(CMAKE_PREFIX_PATH "${pybind11_CMAKE_DIR};${CMAKE_PREFIX_PATH}")

# handle rpath weirdness
if (NOT APPLE)
    set(CMAKE_INSTALL_RPATH $ORIGIN)
endif ()

option(CMAKE_INTERPROCEDURAL_OPTIMIZATION
        "Enable interprocedural optimization for all targets." ON
)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../test/)

# ##############################################################################
# REQUIRED PACKAGES
# ##############################################################################

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

if ("$ENV{SETVARS_COMPLETED}" STREQUAL "1")
    message(
            WARNING
            "Intel OneAPI environment is active. Might lead to issues with correct MPI lib discovery!"
    )
endif ()

find_package (Python3 COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG REQUIRED)

set(HDF5_ENABLE_PARALLEL ON)
find_package(HDF5 REQUIRED COMPONENTS CXX)

find_package(Eigen3 REQUIRED)

find_package(OpenMP REQUIRED)

find_package(MPI REQUIRED)

find_package(FFTW3 REQUIRED COMPONENTS SINGLE DOUBLE LONGDOUBLE MPI)

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
target_compile_definitions(FANS_FANS PUBLIC H5_HAVE_PARALLEL)

if (CMAKE_SYSTEM_PROCESSOR MATCHES "arm|aarch64")
    target_compile_options(FANS_FANS PUBLIC -march=armv8-a+simd+fp+crypto)
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
    target_compile_options(FANS_FANS PUBLIC -mavx2 -mfma)
endif ()

#add_executable(FANS_main)

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
        include/micro.hpp

)

# TODO: when switching to a newer CMake version this will all be done by one
# call to target_link_libraries(FANS_FANS PUBLIC HDF5::HDF5). But CMake 3.16
## does not yet support this.
target_include_directories(FANS_FANS PUBLIC
        ${HDF5_INCLUDE_DIRS}
)

target_include_directories(FANS_FANS PUBLIC ${FFTW3_INCLUDE_DIRS})

# ##############################################################################
# SOURCES
# ##############################################################################

target_sources(FANS_FANS PRIVATE
        src/matmodel.cpp
        src/reader.cpp
)

# target_sources(FANS_main PRIVATE
#         src/micro.cpp
# )

# ##############################################################################
# add pybind moodule
# ##############################################################################
link_directories(/usr/lib/x86_64-linux-gnu/)

include_directories(
    ${pybind11_INCLUDE_DIR}
    ${Python3_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${HDF5_INCLUDE_DIRS}
    ${FFTW3_INCLUDE_DIRS}
)
set(SOURCES src/matmodel.cpp src/reader.cpp)
set(HEADERS include/general.h include/json.hpp include/matmodel.h include/reader.h include/solverCG.h include/solverFP.h include/solver.h include/micro.hpp)
pybind11_add_module(FANS_pybind
${CMAKE_CURRENT_SOURCE_DIR}/src/micro.cpp ${Sources} ${Headers}
)
target_include_directories(FANS_pybind PRIVATE ${Python3_INCLUDE_DIRS})
target_compile_definitions(FANS_pybind PUBLIC H5_HAVE_PARALLEL)
if (CMAKE_SYSTEM_PROCESSOR MATCHES "arm|aarch64")
    target_compile_options(FANS_pybind PUBLIC -march=armv8-a+simd+fp+crypto)
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
    target_compile_options(FANS_pybind PUBLIC -mavx2 -mfma)
endif ()

target_include_directories(FANS_pybind PUBLIC
        ${HDF5_INCLUDE_DIRS}
)

target_include_directories(FANS_pybind PUBLIC ${FFTW3_INCLUDE_DIRS})
#target_include_directories(FANS_pybind PRIVATE ${Python3_INCLUDE_DIRS})
#target_include_directories(FANS_pybind PRIVATE ${pybind11_INCLUDE_DIRS})
target_include_directories(FANS_pybind PRIVATE ${Python3_INCLUDE_DIRS})
target_include_directories(FANS_pybind PRIVATE ${pybind11_INCLUDE_DIRS})

# ##############################################################################
# LINKING
# ##############################################################################

target_link_libraries(FANS_FANS PRIVATE m)

target_link_libraries(FANS_pybind PUBLIC)

# TODO: when switching to a newer CMake version this can all be done by one
# call to target_link_libraries(FANS_FANS PUBLIC HDF5::HDF5). But CMake 3.16
# does not yet support this.
target_link_libraries(FANS_FANS PUBLIC ${HDF5_CXX_LIBRARIES})

target_link_libraries(FANS_FANS PUBLIC MPI::MPI_CXX)

target_link_libraries(FANS_FANS PUBLIC ${FFTW3_LIBRARIES})

target_link_libraries(FANS_FANS PUBLIC Eigen3::Eigen)

target_link_libraries(FANS_pybind PUBLIC ${HDF5_CXX_LIBRARIES})

target_link_libraries(FANS_pybind PUBLIC MPI::MPI_CXX)

target_link_libraries(FANS_pybind PUBLIC ${FFTW3_LIBRARIES})

target_link_libraries(FANS_pybind PUBLIC Eigen3::Eigen)

target_link_libraries(FANS_pybind PRIVATE OpenMP::OpenMP_CXX)

# target_link_libraries(FANS_main PRIVATE FANS::FANS)

# ##############################################################################
# NAMING
# ##############################################################################

set_target_properties(FANS_FANS PROPERTIES OUTPUT_NAME FANS)
set_target_properties(FANS_FANS PROPERTIES EXPORT_NAME FANS)
# set_target_properties(FANS_main PROPERTIES OUTPUT_NAME FANS)

# ##############################################################################
# INSTALLING
# ##############################################################################

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
# install(
#         TARGETS FANS_FANS
#         EXPORT FANS
#         LIBRARY
#         DESTINATION ${CMAKE_INSTALL_LIBDIR}
#         COMPONENT FANS_Runtime
#         NAMELINK_COMPONENT FANS_Development
#         ARCHIVE
#         DESTINATION ${CMAKE_INSTALL_LIBDIR}
#         COMPONENT FANS_Development
#         INCLUDES
#         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FANS
#         PUBLIC_HEADER
#         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FANS
#         COMPONENT FANS_Development
#         PRIVATE_HEADER
#         DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FANS
#         COMPONENT FANS_Development
# )

# install(
#         TARGETS FANS_main
#         RUNTIME
#         DESTINATION ${CMAKE_INSTALL_BINDIR}
#         COMPONENT FANS_Runtime
# )

# install(
#         EXPORT FANS
#         DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FANS
#         NAMESPACE FANS::
#         FILE FANSTargets.cmake
#         COMPONENT FANS_Development
# )

install(
        TARGETS FANS_pybind
        COMPONENT python
        #COMPONENT FANS_Development
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# install(
#         FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/FindFFTW3.cmake
#         DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FANS/modules
#         COMPONENT FANS_Development
# )

# write_basic_package_version_file(
#         ${CMAKE_CURRENT_BINARY_DIR}/FANSConfigVersion.cmake
#         COMPATIBILITY SameMajorVersion
# )
# configure_package_config_file(
#         "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FANSConfig.cmake.in"
#         "${CMAKE_CURRENT_BINARY_DIR}/FANSConfig.cmake"
#         INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FANS"
#         NO_SET_AND_CHECK_MACRO
#         NO_CHECK_REQUIRED_COMPONENTS_MACRO
# )
# install(
#         FILES ${CMAKE_CURRENT_BINARY_DIR}/FANSConfig.cmake
#         ${CMAKE_CURRENT_BINARY_DIR}/FANSConfigVersion.cmake
#         DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FANS
#         COMPONENT FANS_Development
# )

# ##############################################################################
# PACKAGING
# ##############################################################################

if (${PROJECT_IS_TOP_LEVEL})
    add_subdirectory(cmake/packaging)
endif ()
