"""
In this file we initialize the particle system, create coulomb crystal and add new particles into the system
"""

import numpy as np
from parameters import *
from plotting import *
from integration import *
from fileHandling import *

from copy import copy

def RandomVelocity(T=1, mass=ionMass, dim=3):
    maxElement = np.sqrt((Kb*T)/(mass))
    vector = 2 * maxElement * np.random.rand(dim) - maxElement
        
    return vector * (2/f2)#  2 / f2 because of time transformation

def TemperatureToVelocity(T=4, mass=ionMass):
    
    return np.sqrt((8 * Kb * T)/(mass * np.pi)) * (2/f2)

def RandomPosition(maxRadius=0.1 * r0, dim=3):
    
    vector = 2 * maxRadius * np.random.rand(dim) - maxRadius
    
    return vector

def MakeParticleSystem(n=1, m=1, dimension=3, Tion=10, Telectron=0.01):#system with n ions and m electrons
    particles = []

    for i in range(n): 
        position = RandomPosition()
        velocity = RandomVelocity(mass=ionMass, T=Tion) 
        mass = ionMass
        charge = -electronCharge
        
        particle = np.array([position, velocity, mass, charge], dtype=object)
        particles.append(particle)
        
    for i in range(m):
        position = RandomPosition()
        velocity = RandomVelocity(mass=electronMass, T=Telectron)
        mass = electronMass
        charge = electronCharge
        
        particle = np.array([position, velocity, mass, charge], dtype=object)
        particles.append(particle)

    if dimension == 2:
        for particle in particles:                
            particle[0][2] = 0 #setting position in z direction to zero
            particle[1][2] = 0 #setting velocity in z direction to zero
            
    if dimension == 1:
        for particle in particles:                
            particle[0][0] = 0 #setting position in x direction to zero
            particle[1][0] = 0 #setting velocity in x direction to zero
            particle[0][1] = 0 #setting position in y direction to zero
            particle[1][1] = 0 #setting velocity in y direction to zero

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

def MakeCoulombCrystalFromGrid(nCrystal='25', trapParams=np.array([0, 0.4, 0]), tmax=2000*f2/f1 , dt=25):
    
    start = timer()

    n = int(nCrystal)
    #system = MakeParticleSystem(n, 0)
    system = LoadParticleSystem('coulomb_crystals/' + nCrystal) ############temp
    #system = MakeParticleGrid(n=n, T=0.1)

    
    def TotalVelocity(system):
        vTot = 0
        for particle in system:
            v = particle[1]
            vTot = vTot + Norm(v)
            
        return vTot/n
    
    trasholdVel = TemperatureToVelocity(T=0.001, mass=ionMass)   #povodne T=0.01 
    i = -1
    
    nBoost = 10
    while(i < nBoost):
        if len(system) == 0:
            print('NO PARTICLES LEFT')
            return 0
        """
        if (TotalVelocity(system) < trasholdVel) and (i < nBoost - 2):
            for particle in system:
                particle[1] = RandomVelocity(T=5, mass=ionMass)
        """ 
        
        for particle in system:
            if (Norm(particle[1]) < trasholdVel) and (i < nBoost - 3):
                particle[1] = RandomVelocity(T=0.01, mass=ionMass)
        
        rs, vs, stepName, exeTime, energies, system, stability = ODEint(system, 
                                                                        trapParams, tmax, dt, ODESystem=ODESystemEffectiveDamping,
                                                                        Step=StepVerlet, freezeIons=False)
        
        if i < nBoost-3:
            rsFinal = rs
            vsFinal = vs
            exeTimeFinal = 0
            totalEnergyFinal, kineticEnergyFinal, potentialEnergyFinal = [[],[],[]]

        else:
            rsFinal = np.hstack([rsFinal, rs])
            vsFinal = np.hstack([vsFinal, vs])
            totalEnergy, kineticEnergy, potentialEnergy = energies
            totalEnergyFinal = np.hstack([totalEnergyFinal, totalEnergy])
            kineticEnergyFinal = np.hstack([kineticEnergyFinal, kineticEnergy])
            potentialEnergyFinal = np.hstack([potentialEnergyFinal, potentialEnergy])

            
        exeTimeFinal = exeTimeFinal + exeTime       
        i = i + 1
        
    exeTimeFinal = round(exeTimeFinal, 2)

    energiesFinal = [totalEnergyFinal, kineticEnergyFinal, potentialEnergyFinal]
    solution = rsFinal, vsFinal, stepName, exeTimeFinal, energiesFinal, system, 0
    
    SaveParticleSystem(system, 'coulomb_crystals/' + str(n))
    PlotODESolution(solution)
    #PlotFinalPositions(solution)
    PlotEnergy(solution)
    
    stop = timer()
    print(round(stop - start,2),'seconds')

    
      
    
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

def MakeParticle(mass, charge, T=0.01):
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
    
def CreateInitialSystems(mass=electronMass, charge=electronCharge):
    
    systems = []
    rMax = 0.75 * r0
    initRange = 3
    step = rMax / initRange
    v = np.array([1,1,1]) * 5e-7
    vStep = 10
    
    for iX in range(initRange):
        for iZ in range(initRange):
            for iVel in range(initRange):                
                pos = np.array([iX*step, 0, iZ*step])
                vel = v * vStep**iVel
                system = np.array([pos, vel, mass, charge], dtype=object)
                systems.append(np.array([system]))
                
    return np.array(systems)
                