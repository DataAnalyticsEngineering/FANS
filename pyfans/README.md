# pyFANS

Python bindings for FANS, used to drive a large number of FANS micro simulations
from one macro-scale simulation through the
[Micro Manager](https://precice.org/tooling-micro-manager-overview.html) and
preCICE. The class follows the
[Micro Manager API](https://precice.org/tooling-micro-manager-prepare-micro-simulation.html).

pyFANS is a separate Python project that builds against an installed FANS, so
FANS itself has no Python dependency.

## Installing

Install FANS first, then the bindings against it:

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=<prefix>
cmake --build build -j
cmake --install build

CMAKE_PREFIX_PATH=<prefix> pip install ./pyfans
```

Needs Python 3.12 or newer. The result is an `abi3` wheel, so one build works on
every Python from 3.12 onward.

## Usage

```python
from mpi4py import MPI
import numpy as np
import PyFANS

sim = PyFANS.MicroSimulation(sim_id=0, input_file="input.json")
out = sim.solve({"Strains1to3": np.zeros(3), "Strains4to6": np.zeros(3)}, 1.0)
```

The macro strain goes in as Mandel components and comes back as the homogenized
stress and tangent, all as numpy arrays. preCICE vector data is fixed to the mesh
dimension, which is why the six strain components arrive as two 3-vectors and the
tangent comes back as `Cmat1` to `Cmat7`.

See [test_pyfans](../test/test_pyfans/) for a full coupled example.
