
#mpiexec -n 1 ./FANS input_files/test_hyper_elastic.json
#nohup /usr/bin/time -v mpiexec -n 16 ./FANS input_files/3d0_MFL_0_mech.json &
#nohup taskset -c 0 /usr/bin/time -v mpiexec -n 1 ./FANS input_files/berea128_h5_mech.json &
nohup /usr/bin/time -v mpiexec -n 1 ./FANS input_files/sphere_mech.json &
#mpiexec -n 1 ./FANS input_files/berea128_h5_thermal.json