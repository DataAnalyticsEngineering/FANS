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

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if ("$ENV{SETVARS_COMPLETED}" STREQUAL "1")
    message(
            WARNING
            "Intel OneAPI environment is active. Might lead to issues with correct MPI lib discovery!"
    )
endif ()

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

# TODO: move optimization options into presets
target_compile_options(FANS_FANS PUBLIC -O3 -mavx2 -mfma)

add_executable(FANS_main)

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
target_link_libraries(FANS_FANS PUBLIC ${HDF5_CXX_LIBRARIES})

target_link_libraries(FANS_FANS PUBLIC MPI::MPI_CXX)

target_link_libraries(FANS_FANS PUBLIC ${FFTW3_LIBRARIES})

target_link_libraries(FANS_FANS PUBLIC Eigen3::Eigen)

target_link_libraries(FANS_FANS PRIVATE OpenMP::OpenMP_CXX)


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

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
install(
        TARGETS FANS_FANS FANS_main
        EXPORT FANS
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT FANS_Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        COMPONENT FANS_Runtime
        NAMELINK_COMPONENT FANS_Development
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT FANS_Development
        INCLUDES
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FANS COMPONENT FANS_Development
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FANS COMPONENT FANS_Development
        PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FANS COMPONENT FANS_Development
)

install(
        EXPORT FANS
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FANS
        NAMESPACE FANS::
        FILE FANS.cmake
        COMPONENT FANS_Development
)

write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/FANSConfigVersion.cmake
        COMPATIBILITY SameMajorVersion
)
install(FILES FANSConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/FANSConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FANS
        COMPONENT FANS_Development
)

# ##############################################################################
# PACKAGING
# ##############################################################################

if(${PROJECT_IS_TOP_LEVEL})
  add_subdirectory(packaging)
endif()
