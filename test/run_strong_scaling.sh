#!/bin/bash

for ((n=1; n<=16; n++))
do
    echo "Running with $n MPI processes..."
    start_time=$(date +%s.%N)
    mpiexec -n $n ./FANS input_files/sphere256_linear_elastic.json 
    wait
    end_time=$(date +%s.%N)
    runtime=$(echo "$end_time - $start_time" | bc)
    echo "Runtime for $n MPI processes: $runtime seconds"
done
