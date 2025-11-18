#!/bin/bash

# Check if the number of processes is provided as a command line argument
if [ $# -ne 2 ] || [ "$1" != "-n" ]; then
    echo "Usage: $0 -n <num_processes>"
    exit 1
fi

num_processes=$2

# Select executable depending on whether we're inside a conda/pixi env
if [ -n "$CONDA_PREFIX" ]; then
    FANS_EXEC="$CONDA_PREFIX/bin/FANS"
else
    FANS_EXEC="./FANS"
fi

mkdir -p output

FANS_TEST_CASES=(
    # "J2Plasticity"
    # "LinearElastic"
    # "LinearThermal"
    # "PseudoPlastic"
    # "MixedBCs"
    # "CompressibleNeoHookean"
    # "MixedBCs_LargeStrain"
    "problematic"
)
TOTAL_TESTS=${#FANS_TEST_CASES[@]}
passed=0 failed=0 total_time=0 test_num=0

for test_case in "${FANS_TEST_CASES[@]}"; do
    ((test_num++))
    echo "    Start ${test_num}: ${test_case}"
    start_time=$(date +%s.%N)
    if command time -v mpiexec -n $num_processes "$FANS_EXEC" input_files/test_${test_case}.json output/test_${test_case}.h5 > output/test_${test_case}.log 2>&1; then
        status="Passed" && ((passed++))
    else
        status="***Failed" && ((failed++))
    fi
    end_time=$(date +%s.%N)
    elapsed=$(echo "$end_time - $start_time" | bc)
    total_time=$(echo "$total_time + $elapsed" | bc)
    printf "%d/%d Test #%d: %-30s %s %7.2f sec\n" $test_num $TOTAL_TESTS $test_num "$test_case" "$status" $elapsed
done

echo ""
[ $failed -eq 0 ] && echo "100% tests passed, 0 tests failed out of $TOTAL_TESTS" || echo "$(echo "scale=0; 100 * $passed / $TOTAL_TESTS" | bc)% tests passed, $failed tests failed out of $TOTAL_TESTS"
printf "\nTotal Test time (real) = %.2f sec\n" $total_time
exit $failed
