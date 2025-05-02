# Tests

Execute the [`run_tests.sh`](run_tests.sh) file to run tests, which are also examples. For example, to run a linear elastic mechanical homogenization problem for a 6 othonormal load cases on a microstructure image of size `32 x 32 x 32` with a single spherical inclusion,

```bash
mpiexec -n 2 ./FANS input_files/test_LinearElastic.json test_results.h5
```
