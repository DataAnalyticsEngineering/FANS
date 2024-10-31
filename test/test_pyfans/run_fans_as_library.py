"""
Run FANS as a Python callable library.
"""
import PyFANS as fans
import numpy as np

# from mpi4py import MPI

# MPI.Init()

micro = fans.MicroSimulation(1)

macro_data = dict()

macro_data["strains1to3"] = np.array([1, 0, 0])
macro_data["strains4to6"] = np.array([0, 1, 0])

dt = 0.0001

output = micro.solve(macro_data, dt)  # solve FANS for one load vector g0

print(output)

# MPI.Finalize()
