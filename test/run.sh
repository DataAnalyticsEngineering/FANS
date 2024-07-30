nohup /usr/bin/time -v mpiexec --bind-to core -n 12 ./FANS input_files/sphere_ThermalLinear.json test_results.h5 > nohup.log 2>&1 &
