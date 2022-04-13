import numpy as np
from timeit import default_timer as timer
from copy import copy
from intMethods import *
from parameters import *
from coulomb import * 
from equations import *

def DeleteRecombinedParticles(system, rs, vs, rMax):

    for i in range(len(system)):
        deleteThis = []
        if(system[i,3] == 0):
            deleteThis.append(i)
            
    deleteThis = tuple(deleteThis)
    
    system = np.delete(system, deleteThis, 0)
    n = len(system)
    rs = np.delete(rs, deleteThis, 0)
    vs = np.delete(vs, deleteThis, 0)
    rMax = np.delete(rMax, deleteThis)  

    return system, n, rs, vs, rMax               
        
def Recombine(system, i, o):
    system[(i,o),3] = 0   

def ODEint(system, trapParams, tmax=1.3e+2, dt=1e-2, ODESystem=ODESystemExact,  Step=StepRungaKutta):
    
    #dt = GetDt(ODESystem)#get dt depending on system of ODEs you want to solve
    
    particles = copy(system)#don't want to modify initial system
    
    start = timer()#to track real time of the computation
    
    n = len(particles)
    t = 0
    iterations = int((tmax - t) / dt) + 1 #this is number of steps that will by saved in an array
    
    rs = np.zeros((n,iterations,3)) #array of positions in time for all particles
    vs = np.zeros((n,iterations,3)) #array of velocities in time for all particles
    rMax = np.zeros(n) #if rMax[i] > trashold(if it hits the electrode) then we declare particle[i] unstable
    potentialEnergy = np.zeros(iterations) #array of potential energy for the system in time
    kineticEnergy = np.zeros(iterations) #array of kinetic energy for the system in time
    
    stabilityInt = np.zeros(n)#need to modify in case of recombination(in DeleteRecombinedParticles)
        
    for i in range(n): #initial positions and velocities
        rs[i][0] = particles[i][0]
        vs[i][0] = particles[i][1]

    allowRecombination = False
    freezeIons = True        
    
    for k in range(iterations - 1): #loop in time
        
        mass = particles[:,2]
        charge = particles[:,3]
        
        #fCoulomb, potentialCoulomb = CoulombNBody(particles)#test
        
        rMatrix, vMatrix = GetPosVelPairCMS(particles)
        fCoulomb, potentialCoulomb = CoulombNBody(rMatrix, charge)

        aCoulomb = fCoulomb / mass [:,None]
             
        for i in range(n): #loop through particles
            
            r = rs[i][k]
            v = vs[i][k]
            rv = np.array([r,v])
            
            #aCoulomb[i] = np.zeros(3)#test

            
            if(Norm(r) > rMax[i]):
                rMax[i] = Norm(r) 
                
            stabilityInt[i] = stabilityInt[i] + np.dot(r,r)
            
            kineticEnergy[k] = kineticEnergy[k] + 0.5 * mass[i] * Norm(v)**2 * (f2/2)**2
            potentialFromTrap = charge[i] * (V0 + V1 * np.cos(f1 * t) + V2 * np.cos(f2 * t) * (r[0]**2 + r[1]**2 - 2 * r[2]**2) / r0**2)
            potentialEnergy[k] = potentialEnergy[k] + potentialCoulomb[i] + potentialFromTrap 
            
            #potentialEnergy[k] = potentialEnergy[k] - 0.5*const*np.abs(charge[i])*Norm(r)**2 + potentialCoulomb[i] 
                        
            #"""                         
                          
            if allowRecombination:
                
                recombinationEnergy = 1e-12
                recombinationRadius = 1e-10
                
                finerDt = dt
                howMuchFiner = 1           
                
                for o in range(i + 1, n):#another loop thoroug particles -> checking for recombination
                
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
                    rv, t = Step(ODESystem, rv, t, finerDt, aCoulomb[i], mass[i], charge[i], trapParams)
                    
            else:
                if freezeIons:
                    if(charge[i] > 0):
                        t = t + dt
                    else:
                        rv, t = Step(ODESystem, rv, t, dt, aCoulomb[i], mass[i], charge[i], trapParams)
                else:
                    rv, t = Step(ODESystem, rv, t, dt, aCoulomb[i], mass[i], charge[i], trapParams)



            #"""
            #rv, t = Step(ODESystem, rv, t, dt, aCoulomb[i], mass[i], charge[i], trapParams)#test
                
                            
            r, v = rv            
            
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
        #print(rMax[i])
        #if((rMax[i] * 2 / f2) > 0.8 * r0):
        if(rMax[i] > 0.9 * r0):
            stability = stability + 1
            
    stabilityInt = np.sum(stabilityInt) * dt * 2 / (n * tmax * f2)
    #print(stabilityInt)
        
    return rs, vs, Step.__name__, exeTime, energies, particles, stability#stabilityInt
