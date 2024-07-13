nohup /usr/bin/time -v mpiexec --bind-to core -n 16 ./FANS input_files/sphere_HyperElastic.json test_results.h5 > nohup.log 2>&1 &
