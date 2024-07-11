nohup /usr/bin/time -v mpiexec --bind-to core -n 2 ./FANS input_files/test_ThermalLinear.json test_results.h5 > nohup_test_ThermalLinear.log 2>&1 &
nohup /usr/bin/time -v mpiexec --bind-to core -n 2 ./FANS input_files/test_MechLinear.json test_results.h5 > nohup_test_MechLinear.log 2>&1 &
nohup /usr/bin/time -v mpiexec --bind-to core -n 2 ./FANS input_files/test_HyperElastic.json test_results.h5 > nohup_test_HyperElastic.log 2>&1 &
