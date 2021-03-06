"""
In this module handle coulomb interaction between charged particles
"""

import numpy as np
from parameters import *

def GetPosVelPairCMS(particles, freezeIons=False):
    """
    This function takes a system with N particles and returns two NxN matrices. (i,j) element of the first 
    matrix is a distance between i-th and j-th particle in the system. Same for second matrix butwith relative velocities.
    """
    
    n = len(particles)
    
    rs = np.zeros((n,n,3))
    vs = np.zeros((n,n,3))

    #in matrices rs and vs there are stored relative positions and velocities of every aCoulomb[i] of particles
    #we will use these matrices for Coulomb interaction, for recombination and for adaptive time-step
    
    if freezeIons:
        
        electrons = np.array([])
        ions = np.array([])
        
        for i in range(n):
            if particles[i,3] > 0:
                electrons = particles[:i]
                ions = particles[i:]
                break
            
        nElectrons = len(electrons)
        nIons = len(ions)
        
        for i in range(nElectrons):
            r1 = electrons[i,0]
            #v1 = electron[i,1]
            for j in range(i+1, n):
                r2 = particles[j,0] 
                rs[i,j] = r2 - r1
                rs[j,i] = -rs[i,j]
                
                #v2 = 0
                #vs[i,j] = v1 #- v2
                #vs[j,i] = -vs[i,j]
                                
    else:
        for i in range(n):
            r1 = particles[i,0]
            v1 = particles[i,1]
            for j in range(i+1 ,n):
                r2 = particles[j,0]
                rs[i][j] = r2 - r1
                rs[j][i] = -rs[i][j]

                v2 = particles[j,1]
                vs[i][j] = v1 - v2
                vs[j][i] = -rs[i][j]

    return rs, vs
    
def NeedFinerTimeStep(r, v, dt):
    """
    r is vector of relative position between two concrete particles
    v is vector of relative velocity between two concrete particles     
    if projection of (v * timeStep) onto the direction of r is greater than the length of r 
    then we will use finer time step.
    This option is used only when we allow recombinations because it can considerably slow computation.
    """
    a = np.dot(v,r)
    b = np.dot(r,r)

    if(a * dt > 0.7 * b): 
        return True
    else:
        return False     
        

def CoulombForce(r, charges, timeTransformation=True):
    
    """returns Coulomb force and potential energy caused by Coulomb interaction between single pair of particles"""     
    
    if timeTransformation:    
        c = 1 /(4 * np.pi * eps0) * (2/f2)**2 #term (2/f2)**2 is due to time scaling 
    else:
        c = 1 /(4 * np.pi * eps0)
        
    
    force = c * charges * Norm(r)**(-3) * r
    potential = np.dot(force,r)#this is potential energy
    
    return np.array(force), potential

def CoulombNBody(rMatrix, charges, particles, freezeIons=False, timeTransformation=True):
    """
    rMatrix is matrix of relative positions between all pairs of particles
    returns a vector of cumulatice Coulomb force on each particle and a vector of Coulombic potential energy of each particle
    """
    
    n = len(charges)
    forces = np.zeros((n,n,3)) #matrix of Coulomb forces between all pairs of particles
    force = np.zeros((n,3)) #i-th element of this vector is cumulative Coulomb force on i-th particle
    potentials = np.zeros((n,n))
    potential = np.zeros(n)
    

    if freezeIons:
        
        electrons = np.array([])
        ions = np.array([])
        
        for i in range(n):
            if particles[i,3] > 0:
                electrons = particles[:i]
                ions = particles[i:]
                break
            
        nElectrons = len(electrons)
        nIons = len(ions)       
            
        for i in range(nElectrons):
            r1 = electrons[i,0]
            #v1 = electron[i,1]
            for j in range(i+1, n):
                r = rMatrix[i,j]
                chargeProduct = charges[i] * charges[j]
                forces[i][j], potentials[i][j] = CoulombForce(r, chargeProduct, timeTransformation)
                forces[j][i] = -forces[i][j]
                potentials[j][i] = potentials[i][j]
    else:    
        for i in range(n):
            for j in range(i+1 ,n):
                r = rMatrix[i,j]
                chargeProduct = charges[i] * charges[j]
                forces[i][j],potentials[i][j] = CoulombForce(r, chargeProduct, timeTransformation)
                forces[j][i] = -forces[i][j]
                potentials[j][i] = potentials[i][j]          

    for i in range(n):
        for j in range(n):
            force[i] = force[i] + forces[i][j]
            potential[i] = potential[i] + potentials[i][j]

    return force, potential  