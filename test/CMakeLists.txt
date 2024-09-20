set(FANS_TEST_INPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(FANS_TEST_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(FANS_EXECUTABLE $<TARGET_FILE:FANS_main>)

# determine MPI process count. The discretization of the test geometry allows for max. 8 processes.
set(FANS_N_MPI_PROCESSES_MAX 8)
cmake_host_system_information(RESULT FANS_CORES_AVAILABLE QUERY NUMBER_OF_PHYSICAL_CORES)
if (FANS_N_MPI_PROCESSES_MAX LESS FANS_CORES_AVAILABLE)
    set(FANS_N_MPI_PROCESSES ${FANS_N_MPI_PROCESSES_MAX})
else()
    set(FANS_N_MPI_PROCESSES ${FANS_CORES_AVAILABLE})
endif()
message(STATUS "Will use ${FANS_N_MPI_PROCESSES} processes for MPI test cases.")

set(FANS_TEST_CASES
    J2Plasticity
    LinearElastic
    LinearThermal
    PseudoPlastic
)

list(LENGTH FANS_TEST_CASES N_TESTS)
math(EXPR N_TESTS "${N_TESTS} - 1")

foreach(N RANGE ${N_TESTS})
    list(GET FANS_TEST_CASES ${N} FANS_TEST_CASE)

    add_test(
        NAME ${FANS_TEST_CASE}
        COMMAND mpiexec -n ${FANS_N_MPI_PROCESSES} ${FANS_EXECUTABLE} input_files/test_${FANS_TEST_CASE}.json ${FANS_TEST_OUTPUT_DIR}/test_${FANS_TEST_CASE}.h5
        WORKING_DIRECTORY ${FANS_TEST_INPUT_DIR}
    )
endforeach()
