import numpy as np
from parameters import *
from plotting import *
from integration import *
from copy import copy
import csv

def RandomVelocity(T=4, mass=calciumMass, dim=3):
    Kb = 1.3806e-23
    maxElement = np.sqrt((2*Kb*T)/(3*mass))
    vector = 2 * maxElement * np.random.rand(dim) - maxElement
        
    return vector

def RandomPosition(maxRadius=r0, dim=3):
    
    vector = 2 * maxRadius * np.random.rand(dim) - maxRadius
    
    return vector

def MakeParticleSystem(n=1,m=1):#system with n ions and m electrons
    particles = []

    for i in range(n): 
        position = RandomPosition()
        velocity = RandomVelocity() 
        mass = calciumMass
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

    return np.array(particles)

def MakeParticleGrid(n=10):
    
    spacing = 5e-9
    
    pos = np.ones(3) * -(np.cbrt(n)-1) * spacing / 2
    mass = calciumMass
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

def IonCrystal(nIons=10, vMax=2e-1, system=None):
    if (system == None):
        system = MakeParticleGrid(nIons)
    trapParams = np.array([0, 0.45, 0])
        
    tmax = 1e+2
    dt = 1e-1
    
    solution = ODEint(system, trapParams, tmax, dt, ODESystem=ODESystemDampingEff,  Step=StepEulerAdvanced)

    k = 1
    trashold = 10
    
    while((TotalVelocity(system) > vMax) and (k < trashold)):
        solution = ODEint(system, trapParams, tmax, dt, ODESystem=ODESystemDampingEff,  Step=StepEulerAdvanced)
        system = solution[-2]
        k = k + 1
        print(k, TotalVelocity(system))
        
    if(k >= trashold):
        print('unsuccessful')

    #solution = ODEint(system, trapParams, tmax=1e+2, dt=1e-3, ODESystem=ODESystemCrystal,  Step=StepEulerAdvanced)
    #system = solution[-2]  
    PlotODESolution(solution)      
    PlotFinalPositions(solution)              
    
    return system    
    
def AddElectron(system):
    
    r = RandomPosition(maxRadius=1e-9) 
    v = RandomVelocity(1, electronMass)
    m = electronMass
    c = electronCharge
    
    electron = np.array([r, v, m, c], dtype=object)
    newSystem = np.vstack([system,electron])
    
    return newSystem

def AddIon(system):
    
    r = RandomPosition(maxRadius=1e-8) 
    v = RandomVelocity(0.1, calciumMass)
    m = calciumMass
    c = -electronCharge
    
    ion = np.array([r, v, m, c], dtype=object)
    newSystem = np.vstack([system,ion])
    return newSystem

def LoadParticleSystem(fileName='init_system'):
    
    particles = []
    with open(r"data/" + fileName + ".dat") as csvFile:
        csvReader = csv.reader(csvFile, delimiter = "\t")
        
        for row in csvReader:
            
            data = row[0].split(',')
            
            r = np.zeros(3)
            v = np.zeros(3)

            for i in range(3):
                r[i] = float(data[i])                
            for i in range(3):
                v[i] = float(data[i+3]) 
            mass = float(data[6])
            charge = float(data[7])
            
            particle = np.array([r,v,mass,charge], dtype=object)
            particles.append(particle)              
            
    return np.array(particles)

def SaveParticleSystem(system, fileName='init_system'):
    
    n = len(system)
    
    with open(r"data/" + fileName + ".dat","w", newline="") as csvFile:
        csvWriter = csv.writer(csvFile, delimiter = "\t")
        #csvWriter.writerow("Ta toto su nejake data")
        
        for i in range(n):
            row = []
            for j in range(2):
                for k in range(3):
                    row.append(system[i,j][k])
            for j in range(2):
                row.append(system[i,j+2])
            csvWriter.writerow([str(row)[1:-1]])
            
def SaveStabilityDiagram(stability, params):
    
    q1Resol, q2Resol, nParticles = params
    fileName = str(nParticles) + '_particles_' + str(q2Resol) + 'x' + str(q1Resol)
    
    with open(r"data/stability_diagrams/" + fileName + ".dat","w", newline="") as csvFile:
        csvWriter = csv.writer(csvFile, delimiter = "\t")
        #csvWriter.writerow("Ta toto su nejake data")
        for row in stability:
            csvWriter.writerow([str(row)[1:-1]])