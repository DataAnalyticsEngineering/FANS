pybind11_add_module(PyFANS micro.hpp micro.cpp)
target_include_directories(PyFANS PRIVATE ${HDF5_INCLUDE_DIRS} ${FFTW3_INCLUDE_DIRS})
target_link_libraries(PyFANS PRIVATE FANS::FANS)

message(STATUS "PyFANS builds against Python ${Python_VERSION} at ${Python_EXECUTABLE}")

# Not CACHE: a cached value would survive a reconfigure against a different
# Python and install into the old version's site-packages.
if (NOT DEFINED FANS_PYTHON_SITE_DIR)
    set(FANS_PYTHON_SITE_DIR
        "${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages")
endif()

# $ORIGIN/../.. is <prefix>/lib, where libFANS.so lands.
set_target_properties(PyFANS PROPERTIES
    INSTALL_RPATH "$ORIGIN/../.."
    BUILD_WITH_INSTALL_RPATH FALSE)

install(
        TARGETS PyFANS
        LIBRARY
        DESTINATION ${FANS_PYTHON_SITE_DIR}
        COMPONENT FANS_Runtime
)

add_custom_command(
        TARGET PyFANS
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
        $<TARGET_FILE:PyFANS>
        ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_pyfans/$<TARGET_FILE_NAME:PyFANS>
        COMMENT "Create a symlink for FANS python bindings to ${CMAKE_CURRENT_SOURCE_DIR}/../test/"
)
