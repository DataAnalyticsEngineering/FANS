#!/bin/bash

for ((n=1; n<=16; n++))
do
    echo "Running with $n MPI processes..."
    start_time=$(date +%s.%N)
    nohup /usr/bin/time -v mpiexec -n $n ./FANS input_files/sphere256_linear_elastic.json > nohup_$n.out 2>&1 &
    wait
    end_time=$(date +%s.%N)
    runtime=$(echo "$end_time - $start_time" | bc)
    echo "Runtime for $n MPI processes: $runtime seconds"
done
