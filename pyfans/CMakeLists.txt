pybind11_add_module(PyFANS micro.hpp micro.cpp)
target_include_directories(PyFANS PRIVATE ${HDF5_INCLUDE_DIRS} ${FFTW3_INCLUDE_DIRS})
target_link_libraries(PyFANS PRIVATE FANS::FANS)

# Optional integer option (empty by default)
set(MicroSimulationSuffix "" CACHE STRING "Optional numeric variant suffix")

# Validate: must be empty or integer
if(MicroSimulationSuffix STREQUAL "")
    set(CLASS_NAME MicroSimulation)
    set(INPUT_NAME input.json)
    set(MODULE_NAME PyFANS)
else()
    if(NOT MicroSimulationSuffix MATCHES "^[0-9]+$")
        message(FATAL_ERROR "MicroSimulationSuffix must be an integer, got '${MicroSimulationSuffix}'")
    endif()

    set(CLASS_NAME MicroSimulation${MicroSimulationSuffix})
    set(INPUT_NAME input${MicroSimulationSuffix}.json)
    set(MODULE_NAME PyFANS${MicroSimulationSuffix})
    message("Using MicroSimulationSuffix: ${MicroSimulationSuffix}")
endif()

# Add compile definitions
target_compile_definitions(PyFANS PRIVATE
        PyFANS_CLASS_NAME=${CLASS_NAME}
        PyFANS_CLASS_NAME_STR="${CLASS_NAME}"
        PyFANS_INPUT_NAME="${INPUT_NAME}"
        PyFANS_MODULE_NAME=${MODULE_NAME}
)

add_custom_command(
        TARGET PyFANS
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
        $<TARGET_FILE:PyFANS>
        ${CMAKE_CURRENT_SOURCE_DIR}/../test/test_pyfans/$<TARGET_FILE_NAME:PyFANS>
        COMMENT "Create a symlink for FANS python bindings to ${CMAKE_CURRENT_SOURCE_DIR}/../test/"
)
