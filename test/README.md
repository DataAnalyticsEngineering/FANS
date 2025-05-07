# Tests

This directory contains tests for the FANS. Tests serve two purposes:

1. Verify the correct functioning of the code.
2. Provide example problems that demonstrate how FANS can be used.

## Test Structure

The test directory includes:
- `CMakeLists.txt` - Configures tests for CTest integration
- `input_files/` - Example input JSON files for each test case
- Python based validation tests in `pytest/` directory for validating results

## Running Tests

### Using CTest

Tests are configured by the main CMake system. To run all FANS tests:

```bash
cd ../build
ctest fans
```

This will run all test cases and generate result files in `build/test/`.

## Result Validation

After running the tests, verify the results using pytests in the `pytest/` directory with the following command:

```bash
pixi run test
```

Note: The pytest validation expects result files to be in `build/test/` directory, so make sure to run the tests first.

## Available Test Cases

For a 3D microstructure image of resolution `32 x 32 x 32` with a single spherical inclusion, the following test cases are available:
- Linear thermal homogenization problem with isotropic heat conductivity - `test_LinearThermal.json`
- Small strain mechanical homogenization problem with linear elasticity - `test_LinearElastic.json`
- Small strain mechanical homogenization problem with nonlinear pseudoplasticity - `test_PseudoPlastic.json`
- Small strain mechanical homogenization problem with Von-Mises plasticity - `test_J2Plasticity.json`

Each test case has corresponding input JSON files in the `input_files/` directory. Tests can be run individually as example problems. For instance,

```bash
mpiexec -n 2 ./FANS input_files/test_LinearElastic.json test_results.h5
```
