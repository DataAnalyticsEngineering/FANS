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

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpic++)

# with -fno-implicit-templates I get linker errors when using std:: stuff
# TODO: should be developer specific, by using e.g. CMake Presets
set(CMAKE_CXX_FLAGS "-O3 -mavx2 -mfma -fopenmp -lm -DH5_HAVE_PARALLEL -fmax-errors=1")

# handle rpath weirdness
if (NOT APPLE)
    set(CMAKE_INSTALL_RPATH $ORIGIN)
endif ()

option(CMAKE_INTERPROCEDURAL_OPTIMIZATION
        "Enable interprocedural optimization for all targets." ON
)


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY test/)

# ##############################################################################
# REQUIRED PACKAGES
# ##############################################################################

set(HDF5_ENABLE_PARALLEL ON)
find_package(HDF5 REQUIRED COMPONENTS C CXX)

find_package(Eigen3 REQUIRED)

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

# TODO: workaround for not available find_package calls and stuff
# target_include_directories(FANS PRIVATE
#         ${HDF5_INCLUDE_DIRS}
# )
target_include_directories(FANS_FANS PRIVATE
        /usr/local/include
)
target_include_directories(FANS_FANS PRIVATE
        /usr/include
)

target_link_directories(FANS_FANS PRIVATE /usr/lib/x86_64-linux-gnu/)

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

# target_link_libraries(FANS PRIVATE ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES})
target_link_libraries(FANS_FANS PUBLIC HDF5::HDF5)

target_link_libraries(FANS_FANS PUBLIC fftw3 fftw3_mpi)

target_link_libraries(FANS_FANS PUBLIC Eigen3::Eigen)

target_link_libraries(FANS_main PRIVATE FANS::FANS)

# ##############################################################################
# NAMING
# ##############################################################################

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
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FANS
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FANS
        PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/FANS
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
)
