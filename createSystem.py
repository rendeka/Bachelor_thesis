"""
In this file we initialize the particle system, create coulomb crystal and adding new particles into the system
"""

import numpy as np
from parameters import *
from plotting import *
from integration import *
from copy import copy
import csv

def RandomVelocity(T=4, mass=ionMass, dim=3):
    maxElement = np.sqrt((Kb*T)/(mass))
    vector = 2 * maxElement * np.random.rand(dim) - maxElement
        
    return vector * (2/f2)#  2 / f2 because of time transformation

def TemperatureToVelocity(T=4, mass=ionMass):
    
    return np.sqrt((8 * Kb * T)/(mass * np.pi)) * (2/f2)

def RandomPosition(maxRadius=0.4 * r0, dim=3):
    
    vector = 2 * maxRadius * np.random.rand(dim) - maxRadius
    
    return vector

def MakeParticleSystem(n=1,m=1):#system with n ions and m electrons
    particles = []

    for i in range(n): 
        position = RandomPosition()
        velocity = RandomVelocity(mass=ionMass) 
        mass = ionMass
        charge = -electronCharge
        
        particle = np.array([position, velocity, mass, charge], dtype=object)
        particles.append(particle)
        
    for i in range(m):
        position = RandomPosition()
        velocity = RandomVelocity(mass=electronMass)
        mass = electronMass
        charge = electronCharge
        
        particle = np.array([position, velocity, mass, charge], dtype=object)
        particles.append(particle)

    return np.array(particles)

def MakeParticleGrid(n=10):
    
    spacing = 1e-4
    
    pos = np.ones(3) * -(np.cbrt(n)-1) * spacing / 2
    mass = ionMass
    charge = -electronCharge
    
    grid = []
    
    def NextParticle(r, m=mass, c=charge):
        vel = RandomVelocity(1, m) / 5
        pos = copy(r)
        return np.array([pos, vel, m, c], dtype=object)
        
    
    k = 0
    while(k < n):
        grid.append(NextParticle(pos))
        pos[0] = pos[0] + spacing
        k = k + 1
        if(k < n):
            grid.append(NextParticle(pos))
            pos[1] = pos[1] + spacing
            k = k + 1
            if(k < n):
                grid.append(NextParticle(pos))
                pos[2] = pos[2] + spacing
                k = k + 1
                
    return np.array(grid)        
        

def TotalVelocity(system):
    vTot = 0
    for particle in system:
        v = particle[1]
        vTot = vTot + Norm(v)**2
        
    return vTot/len(system)

    
def IonCrystal(nIons=20, vMax=5e-1, system=None, trapParams=np.array([0, 0.4, 0])):
    """making coulomb crystal by molecular dynamics -> solving equations of motion with damping"""
    if (system == None):
        system = MakeParticleGrid(nIons)
            
    vMax = TemperatureToVelocity(T=1e-3, mass=ionMass)#temperature we can get with Doppler cooling is 0.5 miliKelvin for calsium atom
    
    tmax = f2 / f1 
    dt = 1/(100)   
        
    solution = ODEint(system, trapParams, tmax, dt, ODESystem=ODESystemDampingEff,  Step=StepEulerAdvanced)

    k = 1
    trashold = 20
    
    while((TotalVelocity(system) > vMax) and (k < trashold)):
        solution = ODEint(system, trapParams, tmax, dt, ODESystem=ODESystemDampingEff,  Step=StepEulerAdvanced)
        system = solution[-2]
        k = k + 1
        print(k)
        
    if(k >= trashold):
        print('unsuccessful')

    #solution = ODEint(system, trapParams, tmax=1e+2, dt=1e-3, ODESystem=ODESystemCrystal,  Step=StepEulerAdvanced)
    #system = solution[-2]  
    #PlotODESolution(solution)      
    #PlotFinalPositions(solution)              
    
    return system, solution    
    
def AddElectron(system):
    """adds electron to the system"""
    
    r = RandomPosition(maxRadius=1e-4) 
    v = RandomVelocity(10, electronMass)
    m = electronMass
    c = electronCharge
    
    electron = np.array([r, v, m, c], dtype=object)
    newSystem = np.vstack([system,electron])
    
    return newSystem

def AddIon(system):
    """adds ion to the system"""
    
    r = RandomPosition() 
    v = RandomVelocity(10000, ionMass)
    m = ionMass
    c = -electronCharge
    
    ion = np.array([r, v, m, c], dtype=object)
    newSystem = np.vstack([system,ion])
    return newSystem
