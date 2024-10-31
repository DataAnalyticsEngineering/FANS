import PyFANS as fans
import numpy as np

# from mpi4py import MPI

"""
Run FANS from python. For mechanical simulations, import PyFANS. For thermal simulations, use PyFANSTHERMAL.
To be able to run the code, the python bindings must be compiled with the correct flags to initialize MPI.
Move to a build directory and run cmake with the following flags:
cmake -DUSE_MPI=ON ..
"""

# MPI.Init()

micro = fans.MicroSimulation(1)

macro_data = dict()

macro_data["strains1to3"] = np.array([1, 0, 0])
macro_data["strains4to6"] = np.array([0, 1, 0])

dt = 0.0001

output = micro.solve(macro_data, dt)  # solve FANS for one load vector g0

print(output)

# MPI.Finalize()
