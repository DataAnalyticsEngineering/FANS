<<<<<<< HEAD
#!/bin/bash

# Check if the number of processes is provided as a command line argument
if [ $# -ne 2 ] || [ "$1" != "-n" ]; then
    echo "Usage: $0 -n <num_processes>"
    exit 1
fi

num_processes=$2

nohup time -v mpiexec -n $num_processes ./FANS input_files/test_LinearThermalIsotropic.json test_LinearThermalIsotropic.h5 > test_LinearThermalIsotropic.log 2>&1 &

nohup time -v mpiexec -n $num_processes ./FANS input_files/test_LinearElasticIsotropic.json test_LinearElasticIsotropic.h5 > test_LinearElasticIsotropic.log 2>&1 &

nohup time -v mpiexec -n $num_processes ./FANS input_files/test_PseudoPlasticLinearHardening.json test_PseudoPlasticLinearHardening.h5 > test_PseudoPlasticLinearHardening.log 2>&1 &

nohup time -v mpiexec -n $num_processes ./FANS input_files/test_VonMisesPlasticLinearIsotropicHardening.json test_VonMisesPlasticLinearIsotropicHardening.h5 > test_VonMisesPlasticLinearIsotropicHardening.log 2>&1 &
=======
nohup /usr/bin/time -v mpiexec -n 2 ./FANS input_files/test_ThermalLinear.json test_results.h5 > nohup_test_ThermalLinear.log 2>&1 &
nohup /usr/bin/time -v mpiexec -n 2 ./FANS input_files/test_MechLinear.json test_results.h5 > nohup_test_MechLinear.log 2>&1 &
nohup /usr/bin/time -v mpiexec -n 2 ./FANS input_files/test_HyperElastic.json test_results.h5 > nohup_test_HyperElastic.log 2>&1 &
>>>>>>> restructured how to define macrscale loading in the input file
