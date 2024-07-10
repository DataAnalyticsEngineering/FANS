# ##############################################################################
# PACKAGING
# ##############################################################################

set(CPACK_OUTPUT_FILE_PREFIX "${CMAKE_BINARY_DIR}/packages")

set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
set(CPACK_PACKAGE_VENDOR "MIB DAE Stuttgart")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "FANS - Fourier Accelerated Nodal Solver" CACHE STRING "Extended summary.")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.tik.uni-stuttgart.de/DAE/FANS/")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "MIB DAE Stuttgart")

set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_VERBATIM_VARIABLES YES)
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)

# set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_LIST_DIR}/Description.txt)
# set(CPACK_RESOURCE_FILE_WELCOME ${CMAKE_CURRENT_LIST_DIR}/Welcome.txt)
# set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_LIST_DIR}/License.txt)
# set(CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_LIST_DIR}/Readme.txt)

set(CPACK_DEB_COMPONENT_INSTALL ON)
set(CPACK_DEBIAN_ENABLE_COMPONENT_DEPENDS ON)

set(CPACK_DEBIAN_FANS_RUNTIME_PACKAGE_NAME "fans")
set(CPACK_DEBIAN_FANS_DEVELOPMENT_PACKAGE_NAME "fans-dev")

# this option automatically computes the dependencies of shared libraries (by looking at the libs they are themselves
# linked to). Requires 'dpkg-shlibdeps' to be available.
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

# this package is not autodetected but is required for mpi to function properly
set(CPACK_DEBIAN_FANS_RUNTIME_PACKAGE_DEPENDS "openmpi-bin")

# add header packages of dependencies as recommended (they are required to build from the FANS headers).
set(CPACK_DEBIAN_FANS_DEVELOPMENT_PACKAGE_DEPENDS "libhdf5-dev, libopenmpi-dev, libeigen3-dev, libfftw3-dev, libfftw3-mpi-dev")

include(CPack)

cpack_add_component(FANS_Runtime
        DISPLAY_NAME "FANS Runtime"
        DESCRIPTION "FANS shared library and executable"
        REQUIRED
        INSTALL_TYPES Full Developer Minimal
)
cpack_add_component(FANS_Development
        DISPLAY_NAME "FANS Development"
        DESCRIPTION "FANS headers and CMake files"
        DEPENDS FANS_Runtime
        INSTALL_TYPES Full Developer
)
cpack_add_install_type(Full)
cpack_add_install_type(Minimal)
cpack_add_install_type(Developer)