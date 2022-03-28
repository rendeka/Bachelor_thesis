import numpy as np

def Norm(position):
    return np.linalg.norm(position)

electronMass = 9.109e-31
electronCharge = -1.602e-19
mu = 1.66e-27
calciumMass = 40.08 * mu
r0 = 0.5e-3
V0 = 0 #dc potential
V1 = 5#slow ac potential  5
V2 = 100 #fast ac potential 100
f1 = 3e6 * 2 * np.pi
f2 = 2.5e9 * 2 * np.pi

const = 1e-3 #constant for harmonic potential
const2 = 9.5e-1 #constant for damping force