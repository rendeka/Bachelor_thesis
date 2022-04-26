"""
Here's a list of constants (global variables)
"""

import numpy as np

def Norm(position):
    return np.linalg.norm(position)

electronMass = 9.109e-31
electronCharge = -1.602e-19
mu = 1.66e-27
#ionMass = 40.08 * mu#calcium
#ionMass = mu#hydrogen
ionMass = electronMass#positron

r0 = 0.5e-3
V0 = 0 #dc potential
V1 = 5#slow ac potential  5
V2 = 100 #fast ac potential 100

Kb = 1.3806e-23#Boltzman constant

#f1 = 3e6 * 2 * np.pi
#f2 = 2.5e9 * 2 * np.pi

f2 = 6e8 * 2 * np.pi#test
f1 = f2 / 13#test stability diagram

#f2 = 6e8 * 2 * np.pi#test
#f1 = f2 / 170#test trajectory

#f2 = 6e8 * 2 * np.pi#test
#f1 = f2 / 3#test stability diagram



const = 1e-3 #constant for harmonic potential
const2 = 5e-3 #constant for damping force
