#!/bin/bash

for ((n=8; n<=16; n++))
do
    echo "Running with $n MPI processes..."
    start_time=$(date +%s.%N)
    mpiexec --bind-to core -n $n ./FANS input_files/sphere256_MechLinear.json test.h5 > nohup_MechLinear_$n.log 2>&1
    wait
    end_time=$(date +%s.%N)
    runtime=$(echo "$end_time - $start_time" | bc)
    echo "Runtime for $n MPI processes: $runtime seconds"
done
