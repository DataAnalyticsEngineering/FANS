pybind11_add_module(PyFANS micro.hpp micro.cpp)
target_link_libraries(PyFANS PRIVATE FANS::FANS)

pybind11_add_module(PyFANSTHERMAL micro_thermal.hpp micro_thermal.cpp)
target_link_libraries(PyFANSTHERMAL PRIVATE FANS::FANS)

add_custom_command(
        TARGET PyFANS
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
        $<TARGET_FILE:PyFANS>
        ${CMAKE_CURRENT_SOURCE_DIR}/../test/$<TARGET_FILE_NAME:PyFANS>
        COMMENT "Create a symlink for FANS python bindings to ${CMAKE_CURRENT_SOURCE_DIR}/../test/"
)

add_custom_command(
        TARGET PyFANSTHERMAL
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
        $<TARGET_FILE:PyFANSTHERMAL>
        ${CMAKE_CURRENT_SOURCE_DIR}/../test/$<TARGET_FILE_NAME:PyFANSTHERMAL>
        COMMENT "Create a symlink for FANS python bindings to ${CMAKE_CURRENT_SOURCE_DIR}/../test/"
)

option(USE_MPI "Use MPI" OFF)

if (USE_MPI)
        add_definitions(-DUSE_MPI)
endif()