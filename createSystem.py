"""
In this file we initialize the particle system, create coulomb crystal and add new particles into the system
"""

import numpy as np
from parameters import *
from plotting import *
from integration import *
from fileHandling import *

from copy import copy

def RandomVelocity(T=100, mass=ionMass, dim=3):
    maxElement = np.sqrt((Kb*T)/(mass))
    vector = 2 * maxElement * np.random.rand(dim) - maxElement
        
    return vector * (2/f2)#  2 / f2 because of time transformation

def TemperatureToVelocity(T=4, mass=ionMass):
    
    return np.sqrt((8 * Kb * T)/(mass * np.pi)) * (2/f2)

def RandomPosition(maxRadius=0.6 * r0, dim=3):
    
    vector = 2 * maxRadius * np.random.rand(dim) - maxRadius
    
    return vector

def MakeParticleSystem(n=1, m=1, dimension=3):#system with n ions and m electrons
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
        velocity = RandomVelocity(T=4, mass=electronMass)
        mass = electronMass
        charge = electronCharge
        
        particle = np.array([position, velocity, mass, charge], dtype=object)
        particles.append(particle)

    if dimension == 2:
        for particle in particles:                
            particle[0][2] = 0 #setting position in z direction to zero
            particle[1][2] = 0 #setting velocity in z direction to zero

    return np.array(particles)  
def MakeParticleGrid(n=10, T=4):
    
    spacing = 0.5e-4
    
    anchorPos = np.ones(3) * -(np.cbrt(n)-1) * spacing
    mass = ionMass
    charge = -electronCharge
        
    def NextParticle(r, mass=mass, charge=charge):
        vel = RandomVelocity(T=T, mass=mass)
        pos = copy(r)
        return np.array([pos, vel, mass, charge], dtype=object)

    grid = []
    
    iVec = np.array([1,0,0]) * spacing
    jVec = np.array([0,1,0]) * spacing
    kVec = np.array([0,0,1]) * spacing
    
    m = int(np.cbrt(n) + 1) 
    idxMax = 0           
    for i in range(m):
        for j in range(m):
            for k in range(m):
                if idxMax < n:
                    idxMax = idxMax + 1
                    pos = anchorPos + i*iVec + j*jVec + k*kVec
                    grid.append(NextParticle(pos))        
                
    return np.array(grid) 

def MakeCoulombCrystalFromGrid(nCrystal='25', trapParams=np.array([0, 0.4, 0])):
    
    #system = MakeParticleGrid(int(nCrystal))
    system = MakeParticleSystem(int(nCrystal), 0)
    
    dt = 5e-1
    tmax = 5e5 * dt
    
    solution = ODEint(system, trapParams, tmax, dt, ODESystem=ODESystemEffectiveDamping,  Step=StepEulerAdvanced, freezeIons=False)
    
    SaveParticleSystem(system, 'coulomb_crystals/' + str(int(nCrystal)))
    PlotODESolution(solution)
    #PlotFinalPositions(solution)
    #PlotEnergy(solution)
    
      
    
def IonCrystal(nIons=20, vMax=5e-1, system=np.array([]), trapParams=np.array([0, 0.4, 0])):
    """making coulomb crystal by molecular dynamics -> solving equations of motion with damping"""
    if len(system) == 0:
        return system, None
                
    vMax = TemperatureToVelocity(T=1e-3, mass=ionMass)#temperature we can get with Doppler cooling is 0.5 miliKelvin for calsium atom
    
    tmax = 5 * f2 / f1 
    dt = 1/(50)   
        
    solution = ODEint(system, trapParams, tmax, dt, ODESystem=ODESystemEffectiveDamping,  Step=StepEulerAdvanced, freezeIons=False)

    k = 1
    trashold = 20
    
    velocityTest = True    
    while(velocityTest) and (k < trashold):

        solution = ODEint(system, trapParams, tmax, dt, ODESystem=ODESystemEffectiveDamping,  Step=StepEulerAdvanced, freezeIons=False)
        system = solution[-2]
        print(k)
        k = k + 1
        
        velocityTest = False
        for particle in system:
            if Norm(particle[1]) > vMax:
                velocityTest = True

        
    if(k >= trashold):
        print('unsuccessful')

    #solution = ODEint(system, trapParams, tmax=1e+2, dt=1e-3, ODESystem=ODESystemCrystal,  Step=StepEulerAdvanced)
    #system = solution[-2]  
    #PlotODESolution(solution)      
    #PlotFinalPositions(solution)              
    
    return system, solution

def MakeParticle(mass, charge, T=4):
    r = RandomPosition() 
    v = RandomVelocity(mass=mass, T=T)
    
    return np.array([r, v, mass, charge], dtype=object)
    
    
def AddElectron(system, electron=[], T=4):
    """adds electron to the system"""
    if len(electron) == 0:
        r = RandomPosition() 
        v = RandomVelocity(T, electronMass)
        m = electronMass
        c = electronCharge
        
        electron = np.array([r, v, m, c], dtype=object)
    
    if len(system) == 0:
        return np.array([electron], dtype=object)
        
    else:
        return np.vstack([system, electron])
    
def AddIon(system, ion=[], T=4):
    """adds ion to the system"""
    if len(ion) == 0:
        r = RandomPosition() 
        v = RandomVelocity(T, ionMass)
        m = ionMass
        c = -electronCharge
        
        ion = np.array([r, v, m, c], dtype=object)
    
    if len(system) == 0:
        return np.array([ion], dtype=object)
        
    else:
        return np.vstack([system, ion])
    
def MakeCoulombCrystal(nCrystal='20', trapParams=np.array([0, 0.4, 0])):
    """
    CLUSTER
    Makes coulomb crystal using molecular dynamics. Its single argument final number of ions in crystal
    """       
    
    system, solution = IonCrystal(nIons=1, trapParams=trapParams)
    nCrystal = int(nCrystal)
    
    while len(system) < nCrystal:
        
        print('Number of ions: ', len(system))
        system = AddIon(system, T=10)
        system, solution = IonCrystal(system=system, trapParams=trapParams)
        
    system, solution = IonCrystal(system=system, trapParams=trapParams)
    
    SaveParticleSystem(system, 'coulomb_crystals/' + str(int(nCrystal)))
        
    #PlotODESolution(solution)
    #PlotFinalPositions(solution)
    
    print('returning crystal')    
    return system

def TestCoulombCrystal(nCrystal='20', trapParams=np.array([0, 0.4, 0])):
    """Testing whether particles start to move in exact potential"""
    
    system = LoadParticleSystem('coulomb_crystals/' + nCrystal)
    
    for particle in system:
        particle[1] = np.zeros(3)
        
    solution = ODEint(system, trapParams, 1*endTime, timeStep, ODESystem=ODESystemExact,  Step=StepEulerAdvanced, freezeIons=True)      

    print('total velocity of the system', TotalVelocity(system))
    
    SaveParticleSystem(system, 'coulomb_crystals/crystal-evolution_' + nCrystal)