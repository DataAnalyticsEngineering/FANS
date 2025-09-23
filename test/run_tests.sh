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

# Run the jobs serially
command time -v mpiexec -n $num_processes "$FANS_EXEC" input_files/test_LinearThermal.json output/test_LinearThermal.h5 > test_LinearThermal.log 2>&1

command time -v mpiexec -n $num_processes "$FANS_EXEC" input_files/test_LinearElastic.json output/test_LinearElastic.h5 > test_LinearElastic.log 2>&1

command time -v mpiexec -n $num_processes "$FANS_EXEC" input_files/test_PseudoPlastic.json output/test_PseudoPlastic.h5 > test_PseudoPlastic.log 2>&1

command time -v mpiexec -n $num_processes "$FANS_EXEC" input_files/test_J2Plasticity.json output/test_J2Plasticity.h5 > test_J2Plasticity.log 2>&1

command time -v mpiexec -n $num_processes "$FANS_EXEC" input_files/test_MixedBCs.json output/test_MixedBCs.h5 > test_MixedBCs.log 2>&1
