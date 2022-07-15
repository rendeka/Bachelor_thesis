"""
Here's a list of constants (global variables)
"""

import numpy as np

def Norm(position):
    return np.linalg.norm(position)

electronMass = 9.109e-31
electronCharge = -1.602e-19
mu = 1.66e-27
ionMass = 40.08 * mu#calcium
#ionMass = mu#hydrogen
#ionMass = electronMass#positron

r0 = 0.5e-3
V0 = 0 #dc potential
V1 = 5#slow ac potential  5
V2 = 100 #fast ac potential 100

Kb = 1.3806e-23#Boltzman constant
eps0 = 8.8541878e-12
planckConst = 6.62607e-34

#f2 = 2.5e9 * 2 * np.pi
f2 = 1.88e10

f1 = f2 // 833
f2 = f1 * 833

const = 1e-3 #constant for harmonic potential
beta = 1e-5 #constant for damping force

velocityChop = 1100 # trashold value for velocity diagrams

#for a single picture
sizeLabelSmall = 15
sizeTickSmall = 12

# for 2 pictures side by side
sizeLabel = 20
sizeTick = 15

# for 3 pictures side by side
sizeLabelLarge = 30
sizeTickLarge = 20