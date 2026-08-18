"""
Solve one problem directly through the Python bindings
"""

import os
from pathlib import Path
from mpi4py import MPI
import numpy as np
import PyFANS


def main():
    os.chdir(Path(__file__).resolve().parent)

    sim = PyFANS.MicroSimulation(0, False, "input.json", "pyfans-config.json")
    strain = np.array([1.0e-3, 0.0, 0.0, 0.0, 0.0, 0.0])
    out = sim.solve({"Strains1to3": strain[:3], "Strains4to6": strain[3:]}, 1.0)

    stress = np.concatenate([out["Stresses1to3"], out["Stresses4to6"]])
    packed = np.concatenate([out[f"Cmat{i}"] for i in range(1, 8)])
    C = np.zeros((6, 6))
    k = 0
    for i in range(6):
        for j in range(i, 6):
            C[i, j] = C[j, i] = packed[k]
            k += 1

    gap = np.abs(stress - C @ strain).max()
    assert gap < 1.0e-6 * np.abs(stress).max(), (gap, stress, C @ strain)


if __name__ == "__main__":
    main()
