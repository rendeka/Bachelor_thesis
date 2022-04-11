import numpy as np

def Norm(position):
    return np.linalg.norm(position)

electronMass = 9.109e-31
electronCharge = -1.602e-19
mu = 1.66e-27
#ionMass = 40.08 * mu#calcium
ionMass = mu#hydrogen
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



const = 1e-3 #constant for harmonic potential
const2 = 4e-1 #constant for damping force
