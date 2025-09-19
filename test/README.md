# Tests

This directory contains tests for FANS. Tests serve two purposes:

1. Verify the correct functioning of the code.
2. Provide example problems that demonstrate how FANS can be used.

## Test Structure

The test directory includes:

- `CMakeLists.txt` - Configures tests for CTest integration
- `input_files/` - Example input JSON files for each test case
- Python based validation tests in `pytest/` directory for validating results

## Running Tests

### Using CTest

Tests are configured by CMake when FANS is built. If your CMake build folder is named `build`, run the FANS tests in the following way

```bash
cd ../build
ctest fans
```

This will run all test cases and generate result files in `build/test/`.

## Result Validation

After running the tests, the results are verified using pytest. We recommend running pytest via a pre-configured Pixi task,

```bash
pixi run -e dashboard pytest
```

Note: The validation tests expect result files to be in `build/test/` directory, so make sure to run the tests first.

## Available Test Cases

For a 3D microstructure image of resolution `32 x 32 x 32` with a single spherical inclusion, the following test cases are available:

- Linear thermal homogenization problem with isotropic heat conductivity - `test_LinearThermal.json`
- Small strain mechanical homogenization problem with linear elasticity - `test_LinearElastic.json`
- Small strain mechanical homogenization problem with nonlinear pseudoplasticity - `test_PseudoPlastic.json`
- Small strain mechanical homogenization problem with Von-Mises plasticity - `test_J2Plasticity.json`
- Small strain mechanical homogenization problem with linear pseudoplasticity and mixed stress-strain control boundary conditions - `test_MixedBCs.json`

Each test case has corresponding input JSON files in the `input_files/` directory. Tests can be run individually as example problems. For instance,

```bash
mpiexec -n 2 ./FANS input_files/test_LinearElastic.json test_results.h5
```

To quickly visualize the test results, and accompanying XDMF for the HDF5 output can be generated which can be directly opened in ParaView for 3D visualization and analysis:

```bash
pixi run h52xdmf test_results.h5
```
