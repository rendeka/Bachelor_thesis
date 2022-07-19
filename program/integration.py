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

def TotVel(system, nElectrons):
    result = 0
    for i in range(nElectrons):
        result = result + Norm(system[i,1])
    return result
        

def IsExact(ODESystem):
    
    name = ODESystem.__name__.lower()
    
    if 'exact' in name:
        return True
    else:
        return False    

def trapEnergy(ODESystem, trapParams, r, m, t):
    
    a, q1, q2 = trapParams * electronMass / m    
    
    if IsExact(ODESystem):
        return 1/8 * m * f2**2 * (a - 2 * q1 * np.cos(f1 * t) - 2 * q2 * np.cos(f2 * t)) * (r[0]**2 + r[1]**2 - 2 * r[2]**2)

    else:
        if m == electronMass:
            return 1/8 * m * (a + q2**2 / 2) * f2**2 * (r[0]**2 + r[1]**2 + 2 * r[2]**2)
        else:
            return 1/8 * m * (a + q1**2 / 2 * (f2/f1)**2 + q2**2 / 2) * f2**2 * (r[0]**2 + r[1]**2 + 2 * r[2]**2)
            
def RearrangeSystem(particles):
    """
    we want to rearrange list of particles so we can hadle them more effectively
    """ 
    n = len(particles)
    nElectrons = 0
    
    electrons = []
    ions = []    
    
    for particle in particles:
        if particle[3] < 0:
            electrons.append(particle)
            nElectrons = nElectrons + 1
        else:
            ions.append(particle)
            
    nIons = n - nElectrons
    
    return np.array(electrons + ions), n, nElectrons, nIons

def DeleteRecombinedParticles(system, rs, vs, rMax):
    """removes recombined particles from the system"""
    
    deleteThis = []
    for i in range(len(system)):
        if(system[i,3] == 0): #if charge of a particle is zero than remoze it from the system
            deleteThis.append(i)
            
    deleteThis = tuple(deleteThis)
    
    system = np.delete(system, deleteThis, 0)
    n = len(system)
    #rs = np.delete(rs, deleteThis, 0)
    #vs = np.delete(vs, deleteThis, 0)
    #rMax = np.delete(rMax, deleteThis)  

    return system, n, rs, vs, rMax               
        
def Recombine(system, i, o):
    """marking recombined particles by setting their charge to zero. Later these particles will be removed from the system"""
    system[(i,o),3] = 0   

def ODEint(system, trapParams, tmax=1.3e+2, dt=1e-2, ODESystem=ODESystemExact,
           Step=StepEulerAdvanced, freezeIons=False, velocityDiagram=False, kSecular=20):
    """
    integrating equations of motion
    """
    
    particles = copy(system)#don't want to modify initial system
    particles, n, nElectrons, nIons = RearrangeSystem(particles)
    
    if freezeIons: #with freezed ions we will loop only through electrons
        n = nElectrons
        
    a, q1, q2 = trapParams
    #dt = GetDt(ODESystem)#get dt depending on system of ODEs you want to solve
    
    if (nIons == 0) or freezeIons:
        if q2 > 0:
            if IsExact(ODESystem):    
                tmax = kSecular * np.sqrt(2) / q2
            else:
                tmax = kSecular * np.sqrt(8) / (f2 * q2)
                
    else:
        if q1 > 0:
            if IsExact(ODESystem):    
                tmax = kSecular * np.sqrt(2) / q1
            else:
                dt = dt * ionMass / electronMass
                tmax = 500*kSecular * np.sqrt(2) / q1

                #tmax = kSecular * np.sqrt(8) / (f1 * q1)
                print('what are you doing?')
                
                
    t = np.zeros(n)
    iterations = int(tmax / dt) + 1 #this is number of steps that will by saved in an array
    
    
    rs = np.zeros((n,iterations,3)) #array of positions in time for all particles
    vs = np.zeros((n,iterations,3)) #array of velocities in time for all particles
    rMax = np.zeros(n) #if rMax[i] > trashold(if it hits the electrode) then we declare particle[i] unstable
    potentialEnergy = np.zeros(iterations) #array of potential energy for the system in time
    kineticEnergy = np.zeros(iterations) #array of kinetic energy for the system in time
            
    for i in range(n): #initial positions and velocities
        rs[i][0] = particles[i][0]
        vs[i][0] = particles[i][1]

    allowRecombination = False
    timeTransform = IsExact(ODESystem)
    
    velInit = TotVel(particles, nElectrons)
    velFinal = 0
    
    start = timer()#to track real time of the computation

    for k in range(iterations - 1): #loop in time
        
        mass = particles[:,2]
        charge = particles[:,3]
                
        rMatrix, vMatrix = GetPosVelPairCMS(particles, freezeIons)
        fCoulomb, potentialCoulomb = CoulombNBody(rMatrix, charge, particles, freezeIons, timeTransform)

        aCoulomb = fCoulomb / mass [:,None]
             
        for i in range(n): #loop through particles
            
            r = rs[i][k]
            v = vs[i][k]
            rv = np.array([r,v])
            
            if(Norm(r) > rMax[i]): #tracking the most distant point in trajectory
                rMax[i] = Norm(r)
                #if rMax[i] > r0: #temporary for making coulomb crystal
                #    particles[i,3] = 0
                            
            kineticEnergy[k] = kineticEnergy[k] + 0.5 * mass[i] * Norm(v)**2 * (f2/2)**2
            potentialFromTrap = trapEnergy(ODESystem, trapParams, r, mass[i], t[i]) 
            potentialEnergy[k] = potentialEnergy[k] + potentialFromTrap + potentialCoulomb[i]
            
            
            potentialEnergy[k] = potentialEnergy[k] + potentialCoulomb[i]
            kineticEnergy[k] = 0
                                                              
            if allowRecombination:
                
                recombinationEnergy = 1e-12
                recombinationRadius = 1e-10
                
                finerDt = dt
                howMuchFiner = 1           
                
                for o in range(i + 1, n):#another loop throughout particles -> checking for recombination
                
                    finerDt = dt
                    howMuchFiner = 1#now for skipping data saving
                
                    rCMS = rMatrix[i,o] #CMS -> central mass system(different for every pair of particles)
                    vCMS = vMatrix[i,0] + (aCoulomb[i] - aCoulomb[o]) * dt
                    massCMS = (mass[i] * mass[o]) / (mass[i] + mass[o])
                    #"""
                    if NeedFinerTimeStep(rCMS, vCMS, dt):
                        howMuchFiner = 100
                        finerDt = dt / howMuchFiner
                    #"""
                    CMSKineticEnergy = 0.5 * massCMS * Norm(vCMS)**2 * np.sign(charge[i] * charge[o])
                    
                    #if((Norm(rCMS) < recombinationRadius) and (CMSKineticEnergy < recombinationEnergy)):
                    #    particles = Recombine(particles,i,o)
                        
                    
                for _ in range(howMuchFiner):
                    rv, t[i] = Step(ODESystem, rv, t[i], finerDt, aCoulomb[i], mass[i], charge[i], trapParams)
                    
            else:
                rv, t[i] = Step(ODESystem, rv, t[i], dt, aCoulomb[i], mass[i], charge[i], trapParams)
                           
            r, v = rv
            
            if np.isposinf(np.dot(r,r)) or np.isposinf(np.dot(v,v)): #if position or velocity is too large we stop integrating
                if velocityDiagram:
                    stability = velocityChop
                else:
                    stability = nElectrons
                end = timer()#to track real time of the computation
                exeTime = round((end - start), 2)
                energy = kineticEnergy + potentialEnergy
                energies = [energy, kineticEnergy, potentialEnergy]
                
                return rs, vs, Step.__name__, exeTime, energies, particles, stability
                 
            rs[i][k+1] = r
            vs[i][k+1] = v
            
        for i in range(n):
            particles[i][0] = rs[i][k+1]
            particles[i][1] = vs[i][k+1]
            
        if allowRecombination:
            particles, n, rs, vs, rMax = DeleteRecombinedParticles(particles, rs, vs, rMax)
            if len(particles) == 0:
                return rs, vs, Step.__name__, timer()-start, [[],[],[]], particles, 1000
            
        velFinal = velFinal + TotVel(particles, nElectrons)

    
    end = timer()#to track real time of the computation
    exeTime = round((end - start), 2)
    
    kineticEnergy = kineticEnergy[:-1]#just deleting last element
    potentialEnergy = potentialEnergy[:-1]    
    
    energy = kineticEnergy + potentialEnergy
    energies = [energy, kineticEnergy, potentialEnergy]
    
    velFinal = velFinal / (iterations - 1)
    
    if velocityDiagram:
        stability = round((velFinal / velInit), 5)
        if stability > velocityChop-100:
            stability = velocityChop
    else:
        stability = 0
        
        for i in range(n):
            if(rMax[i] > 0.8 * r0) and (charge[i] < 0): # condition (charge[i] < 0) is here for the case of freezed ions
                stability = stability + 1
        
        
    if freezeIons:
        for particle in particles:
            if particle[3] > 0:
                particle[3] = 0
        particles, n, rs, vs, rMax = DeleteRecombinedParticles(particles, rs, vs, rMax)      
    
    """WARNING! other parts of the program expects that the last value that this function (ODEint) returns is the stability parameter"""    
    return rs, vs, Step.__name__, exeTime, energies, particles, stability

"""
def PlotPotential(ODESystem=ODESystemExact, trapParams=np.array([0, 0.14, 0.37]), m=electronMass, t=2.2):
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable

    
    x_vals = np.linspace(-2, 2, 400)
    z_vals = np.linspace(0, 2, 400)
    X, Z = np.meshgrid(x_vals, z_vals)    
            
    def Phi(X,Z):        
        return trapEnergy(ODESystem, trapParams, np.array([X,0,Z], dtype=object), m, t)
    
    
    fig = plt.figure()        
    plt.contourf(X, Z, Phi(X, Z))
"""    
    
"""
    fig, ax = plt.subplots()
    
    vmin, vmax = 0, 1
    levels = 5
    level_boundaries = np.linspace(vmin, vmax, levels + 1)
    
    quadcontourset = ax.contourf(
        X, Z, Phi(X,Z),
        levels,
        vmin=vmin, vmax=vmax
    )
    
    
    fig.colorbar(
        ScalarMappable(norm=quadcontourset.norm, cmap=quadcontourset.cmap),
        ticks=range(vmin, vmax+5, 5),
        boundaries=level_boundaries,
        values=(level_boundaries[:-1] + level_boundaries[1:]) / 2,
        extend='max',
    )

    
    plt.show()
"""