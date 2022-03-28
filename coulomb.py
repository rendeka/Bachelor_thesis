import numpy as np

def Norm(position):
    return np.linalg.norm(position)

#CoulombForce -> Coulomb interaction between two particles
#"""
def CoulombForce(particle1, particle2):
    
    r1 = particle1[0]
    r2 = particle2[0]
    r = r1 - r2
    
    charges = particle1[3] * particle2[3]
    
    c = 4 * np.pi * 8.854e-12
    force = 1/c * charges * Norm(r)**(-3) * r
    potential = np.dot(force,r)#this is potential energy
    
    return np.array(force), potential

#CoulombNBody -> Coulomb interaction between all particles + potential energy
def CoulombNBody(particles):
    n = len(particles)
    forces = np.zeros((n,n,3))
    force = np.zeros((n,3))
    potentials = np.zeros((n,n))
    potential = np.zeros(n)
    
    for i in range(n):
        p1 = particles[i]
        for j in range(i+1 ,n):
            p2 = particles[j]
            forces[i][j],potentials[i][j] = CoulombForce(p1, p2)
            forces[j][i] = -forces[i][j]
            potentials[j][i] = potentials[i][j]
            

    for i in range(n):
        for j in range(n):
            force[i] = force[i] + forces[i][j]
            potential[i] = potential[i] + potentials[i][j]

    return np.array(force), potential    
#"""

"""

def GetPosVelPairCMS(particles):
    
    n = len(particles)
    
    rs = np.zeros((n,n,3))
    vs = np.zeros((n,n,3))

    #in matrices rs and vs there are stored relative positions and velocities of every aCoulomb[i] of particles
    #we will use these matrices for Coulomb interaction, for recombination and for adaptive time-step
    
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
    a = np.dot(v,r)
    b = np.dot(r,r)
    
    if(a * dt > 0.6 * b): 
        return True
    else:
        return False     
        

def CoulombForce(r, charges):
    
    c = 4 * np.pi * 8.854e-12
    force = 1/c * charges * Norm(r)**(-3) * r
    potential = np.dot(force,r)#this is potential energy
    
    return np.array(force), potential

#CoulombNBody -> Coulomb interaction between all particles + potential energy
def CoulombNBody(rMatrix, charges):
    n = len(charges)
    forces = np.zeros((n,n,3))
    force = np.zeros((n,3))
    potentials = np.zeros((n,n))
    potential = np.zeros(n)
    
    for i in range(n):
        for j in range(i+1 ,n):
            r = rMatrix[i,j]
            chargeProduct = charges[i] * charges[j]
            forces[i][j],potentials[i][j] = CoulombForce(r, chargeProduct)
            forces[j][i] = -forces[i][j]
            potentials[j][i] = potentials[i][j]
            

    for i in range(n):
        for j in range(n):
            force[i] = force[i] + forces[i][j]
            potential[i] = potential[i] + potentials[i][j]

    return force, potential  
"""   