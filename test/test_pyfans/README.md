# Test pyFANS

Test pyFANS as standalone library called from a Python script.

## Build pyFANS

Configure the FANS CMake build with the variable `FANS_LIB` set to `ON`.

## Run the test

Run

```bash
python3 run_fans_as_library.py
```

The script creates a pyFANS object and calls the `solve()` method. The script only checks if the pyFANS object is created and the solve function is callable. The result is not checked for correctness.
