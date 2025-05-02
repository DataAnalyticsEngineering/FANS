# Tests

**TODO: Add a few sentences on what is being tested** Tests are run by executing the [`run_tests.sh`](run_tests.sh) scripts. The tests can also be run individually, as they are basically example problems. For example, to run a linear elastic mechanical homogenization problem for 6 othonormal load cases on a microstructure image of size `32 x 32 x 32` with a single spherical inclusion, in parallel, run

```bash
mpiexec -n 2 ./FANS input_files/test_LinearElastic.json test_results.h5
```
