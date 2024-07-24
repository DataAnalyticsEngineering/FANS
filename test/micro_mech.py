import PyFANS
import numpy as np

micro = PyFANS.MicroSimulation(1) # initialize the micro simulation with id 1

micro.initialize() # initialize FANS

g0 = {"g0": np.array([1,0,0,0,0,0])}

dt = 0.0001

output = micro.solve(g0, dt) # solve FANS for one load vector g0

print(output)
    