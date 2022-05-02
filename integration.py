"""
Main task of this module is to integrate given equations of motion for system of charged particles.
We can also track energy bilance of the system, allow or disallow recombinations and freeze heavy ions
for solving equations of motion for electrons in a bigger Coulomb crystal
"""

import numpy as np
from timeit import default_timer as timer
from copy import copy
from intMethods import *
from parameters import *
from coulomb import * 
from equations import *

def DeleteRecombinedParticles(system, rs, vs, rMax):
    """removes recombined particles from the system"""

    for i in range(len(system)):
        deleteThis = []
        if(system[i,3] == 0): #if charge of a particle is zero than remoze it from the system
            deleteThis.append(i)
            
    deleteThis = tuple(deleteThis)
    
    system = np.delete(system, deleteThis, 0)
    n = len(system)
    rs = np.delete(rs, deleteThis, 0)
    vs = np.delete(vs, deleteThis, 0)
    rMax = np.delete(rMax, deleteThis)  

    return system, n, rs, vs, rMax               
        
def Recombine(system, i, o):
    """marking recombined particles by setting their charge to zero. Later these particles will be removed from the system"""
    system[(i,o),3] = 0   

def ODEint(system, trapParams, tmax=1.3e+2, dt=1e-2, ODESystem=ODESystemExact,  Step=StepRungaKutta, freezeIons=False):
    """
    integrating equations of motion
    ------
    """
    
    #dt = GetDt(ODESystem)#get dt depending on system of ODEs you want to solve
    
    particles = copy(system)#don't want to modify initial system
        
    start = timer()#to track real time of the computation
    
    n = len(particles)
    nIons = 0
    for particle in particles:
        if particle[3] > 0:
            nIons = nIons + 1
    nElectrons = n - nIons
    
    t = np.zeros(n)
    iterations = int((tmax - t[0]) / dt) + 1 #this is number of steps that will by saved in an array
    
    rs = np.zeros((n,iterations,3)) #array of positions in time for all particles
    vs = np.zeros((n,iterations,3)) #array of velocities in time for all particles
    rMax = np.zeros(n) #if rMax[i] > trashold(if it hits the electrode) then we declare particle[i] unstable
    potentialEnergy = np.zeros(iterations) #array of potential energy for the system in time
    kineticEnergy = np.zeros(iterations) #array of kinetic energy for the system in time
            
    for i in range(n): #initial positions and velocities
        rs[i][0] = particles[i][0]
        vs[i][0] = particles[i][1]

    allowRecombination = False
    
    for k in range(iterations - 1): #loop in time
        
        mass = particles[:,2]
        charge = particles[:,3]
                
        rMatrix, vMatrix = GetPosVelPairCMS(particles)
        fCoulomb, potentialCoulomb = CoulombNBody(rMatrix, charge)

        aCoulomb = fCoulomb / mass [:,None]
             
        for i in range(n): #loop through particles
            
            r = rs[i][k]
            v = vs[i][k]
            rv = np.array([r,v])
            
            if(Norm(r) > rMax[i]): #tracking the most distant point in trajectory
                rMax[i] = Norm(r)
                            
            """
            kineticEnergy[k] = kineticEnergy[k] + 0.5 * mass[i] * Norm(v)**2 * (f2/2)**2
            potentialFromTrap = mass[i] * f2**2 * r0**2 * 1/4 * (trapParams[0] / 2 - trapParams[1] * np.cos(f1 * t[i]) - trapParams[2] * np.cos(f2 * t[i]) * (r[0]**2 + r[1]**2 - 2 * r[2]**2))
            potentialEnergy[k] = potentialEnergy[k] + potentialCoulomb[i] + potentialFromTrap
            """
                                                              
            if allowRecombination:
                
                recombinationEnergy = 1e-12
                recombinationRadius = 1e-10
                
                finerDt = dt
                howMuchFiner = 1           
                
                for o in range(i + 1, n):#another loop throughout particles -> checking for recombination
                
                    finerDt = dt
                    howMuchFiner = 1
                
                    rCMS = rMatrix[i,o] #CMS -> central mass system(different for every pair of particles)
                    vCMS = vMatrix[i,0] + (aCoulomb[i] - aCoulomb[o]) * dt
                    massCMS = (mass[i] * mass[o]) / (mass[i] + mass[o])

                    if NeedFinerTimeStep(rCMS, vCMS, dt):
                        howMuchFiner = 100
                        finerDt = dt / howMuchFiner
                    
                    if(charge[i] * charge[o] < 0):
                        CMSKineticEnergy = 0.5 * CMSmass * Norm(CMSv)**2
                        
                        if((rCMS < recombinationRadius) and (CMSKineticEnergy < recombinationEnergy)):
                            particles = Recombine(particles,i,o)
                    
                for _ in range(howMuchFiner):
                    rv, t[i] = Step(ODESystem, rv, t[i], finerDt, aCoulomb[i], mass[i], charge[i], trapParams)
                    
            else:
                if freezeIons:
                    if(charge[i] > 0):
                        t[i] = t[i] + dt
                    else:
                        rv, t[i] = Step(ODESystem, rv, t[i], dt, aCoulomb[i], mass[i], charge[i], trapParams)
                else:
                    rv, t[i] = Step(ODESystem, rv, t[i], dt, aCoulomb[i], mass[i], charge[i], trapParams)
                           
            r, v = rv
            
            if np.isposinf(np.dot(r,r)) or np.isposinf(np.dot(v,v)): #if position or velocity is too large we stop integrating
                stability = nElectrons
                return rs, None, Step.__name__, None, None, particles, stability
                 
            rs[i][k+1] = r
            vs[i][k+1] = v
            
        for i in range(n):
            particles[i][0] = rs[i][k+1]
            particles[i][1] = vs[i][k+1]
            
        if allowRecombination:
            particles, n, rs, vs, rMax = DeleteRecombinedParticles(particles, rs, vs, rMax)
    
    end = timer()#to track real time of the computation
    exeTime = round((end - start), 2)
    
    kineticEnergy = kineticEnergy[:-1]#just deleting last element
    potentialEnergy = potentialEnergy[:-1]    
    
    energy = kineticEnergy + potentialEnergy
    energies = [energy, kineticEnergy, potentialEnergy]
    
    stability = 0
    for i in range(n):
        if(rMax[i] > 0.8 * r0) and (charge[i] < 0): # condition (charge[i] < 0) is here for the case of freezed ions
            stability = stability + 1

    
    """WARNING! other parts of the program expects that the last value that this function (ODEint) returns is the stability parameter"""    
    return rs, vs, Step.__name__, exeTime, energies, particles, stability

