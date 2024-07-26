import PyFANS as fans # alternatively PyFANSTHERMAL
import numpy as np

"""
Run FANS from python. For mechanical simulations, import PyFANS. For thermal simulations, use PyFANSTHERMAL.
To be able to run the code, the python bindings must be compiled with the correct flags to initialize MPI.
Move to a build directory and run cmake with the following flags:
cmake -DUSE_MPI=ON ..
"""

micro = fans.MicroSimulation(1)


g0 = {"g0": np.array([1,0,0,0,0,0])}

dt = 0.0001

output = micro.solve(g0, dt) # solve FANS for one load vector g0

print(output)
