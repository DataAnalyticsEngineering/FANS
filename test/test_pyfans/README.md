# Test pyFANS

Test pyFANS as standalone library called from a Python script.

## Build pyFANS

Configure the FANS CMake build with the variable `FANS_LIB` set to `ON`.

## Dependencies

Install the following dependencies

## Run the test

The test runs a dummy macro problem (unit cube) which is coupled via preCICE to the Micro Manager. The Micro Manager controls micro simulations created using pyFANS. Run the test by running

```bash
python macro-cube.py & micro-manager-precice micro-manager-config.json
```
