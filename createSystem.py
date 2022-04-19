import numpy as np
from parameters import *
from plotting import *
from integration import *
from copy import copy
import csv

def RandomVelocity(T=4, mass=ionMass, dim=3):
    maxElement = np.sqrt((Kb*T)/(mass)) * (2/f2)#the term * (2/f2) is due to timr scaling
    vector = 2 * maxElement * np.random.rand(dim) - maxElement
        
    return vector# * f2 / 2# * f2 / 2 because of time transformation

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
    
    spacing = 5e-5
    
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

    
def IonCrystal(nIons=10, vMax=5e-1, system=None):
    if (system == None):
        system = MakeParticleGrid(nIons)
    trapParams = np.array([0, 0.4, 0])
    
    vMax = TemperatureToVelocity(T=1e-3, mass=ionMass)#temperature we can get with Doppler cooling is 0.5 miliKelvin for calsium atom
    
    tmax = f2 / f1 #in normal time scale 5 * (1 / f1) -> five slower periods
    dt = 1/(100)   
        
    solution = ODEint(system, trapParams, tmax, dt, ODESystem=ODESystemDampingEff,  Step=StepEulerAdvanced)

    k = 1
    trashold = 25
    
    while((TotalVelocity(system) > vMax) and (k < trashold)):
        solution = ODEint(system, trapParams, tmax, dt, ODESystem=ODESystemDampingEff,  Step=StepEulerAdvanced)
        system = solution[-2]
        k = k + 1
        print(k, TotalVelocity(system))
        
    if(k >= trashold):
        print('unsuccessful')

    #solution = ODEint(system, trapParams, tmax=1e+2, dt=1e-3, ODESystem=ODESystemCrystal,  Step=StepEulerAdvanced)
    #system = solution[-2]  
    #PlotODESolution(solution)      
    #PlotFinalPositions(solution)              
    
    return system, solution    
    
def AddElectron(system):
    
    r = RandomPosition(maxRadius=1e-9) 
    v = RandomVelocity(1, electronMass)
    m = electronMass
    c = electronCharge
    
    electron = np.array([r, v, m, c], dtype=object)
    newSystem = np.vstack([system,electron])
    
    return newSystem

def AddIon(system):
    
    r = RandomPosition() 
    v = RandomVelocity(0.1, ionMass)
    m = ionMass
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
    
    q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, time, f1, f2 = params
    nIons, nElectrons = nParticles
    fileName = str(int(nIons)) + '_ions_' + str(int(nElectrons)) + '_electrons_' + 'q1_' + str(q1Start) + '-' + str(q1Stop) + '_q2_' + str(q2Start) + '-' + str(q2Stop)+ '_' + str(int(q2Resol)) + 'x' + str(int(q1Resol)) + '_' + str(int(f2/f1))
    
    with open(r"data/stability_diagrams/" + fileName + ".dat","w", newline="") as csvFile:
        csvWriter = csv.writer(csvFile, delimiter = "\t")
        #csvWriter.writerow("Ta toto su nejake data")
        for row in stability:
            csvWriter.writerow([str(row)[1:-1]])
            
def ParseFileName(fileName='0_ions_1_electrons_q1_0-0.06_q2_0-0.48_700x700_13'):
    parseFileName = fileName.split('_')
    
    nIons = int(parseFileName[0])
    nElectrons = int(parseFileName[2])
    
    parseQ1 = parseFileName[5].split('-')
    q1Start = float(parseQ1[0])
    q1Stop = float(parseQ1[1])
    
    parseQ2 = parseFileName[7].split('-')
    q2Start = float(parseQ2[0])
    q2Stop = float(parseQ2[1])
    
    parseResol = parseFileName[8].split('x')
    q1Resol = int(parseResol[1])
    q2Resol = int(parseResol[0])
    
    eta = parseFileName[9]
    
    return np.array([q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nIons, nElectrons, eta],dtype=object)
            
def LoadStabilityDiagram(fileName='0_ions_1_electrons_q1_0-0.06_q2_0-0.48_700x700_13'):
    
    q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nIons, nElectrons, eta = ParseFileName(fileName)
    
    stability = np.zeros((q1Resol,q2Resol))
    
    with open(r"data/stability_diagrams/" + fileName + ".dat") as csvFile:
        csvReader = csv.reader(csvFile, delimiter = "\t")
        
        i = 0        
        for row in csvReader:
            
            dataRaw = row[0].split(' ')
            data = []
            for k in range(len(dataRaw)):
                prepareDataValue = dataRaw[k].split('.')
                data.append(int(prepareDataValue[0]))

            for j in range(len(data)):
                stability[i,j] = data[j]
                
            i = i + 1
            
    params = np.array([q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, np.array([nIons, nElectrons]), eta, None, None], dtype=object)

    return stability, params
    
def SaveTriangles(triangleUnstable, triangleStable, params):
    
    q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, eta, f1, f2 = params
    nIons, nElectrons = nParticles
    n = nIons + nElectrons
    
    nUnstable = len(triangleUnstable) // 3
    unstables = []
    for i in range(nUnstable):
        unstables.append([triangleUnstable[i*3:i*3 + 3], n])
        
    
    nStable = len(triangleStable) // 3
    stables = []
    for i in range(nStable):
        stables.append([triangleStable[i*3:i*3 + 3], n])
    
    fileNameStable = 'stable_' + str(int(nIons)) + '_ions_' + str(int(nElectrons)) + '_electrons_' + 'q1_' + str(q1Start) + '-' + str(q1Stop) + '_q2_' + str(q2Start) + '-' + str(q2Stop) + '_' + str(int(eta))
    fileNameUnstable = 'unstable_' + str(int(nIons)) + '_ions_' + str(int(nElectrons)) + '_electrons_' + 'q1_' + str(q1Start) + '-' + str(q1Stop) + '_q2_' + str(q2Start) + '-' + str(q2Stop) + '_' + str(int(eta))
    
    with open(r"data/triangles/" + fileNameUnstable + ".dat","w", newline="") as csvFile:
        csvWriter = csv.writer(csvFile, delimiter = "\t")
        #csvWriter.writerow("Ta toto su nejake data")          
        for unstable in unstables:

            unstable = str(unstable)
            unstable = unstable.replace('[','')
            unstable = unstable.replace(']','')
            csvWriter.writerow([unstable])
            
    with open(r"data/triangles/" + fileNameStable + ".dat","w", newline="") as csvFile:
        csvWriter = csv.writer(csvFile, delimiter = "\t")
        #csvWriter.writerow("Ta toto su nejake data")
        for stable in stables:

            stable = str(stable)
            stable = stable.replace('[','')
            stable = stable.replace(']','')
            csvWriter.writerow([stable])
            
            
def LoadTriangles(params):
    
    q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, eta, f1, f2 = params
    nIons, nElectrons = nParticles
    n = nIons + nElectrons
    
    stableFileName = 'stable_' + str(int(nIons)) + '_ions_' + str(int(nElectrons)) + '_electrons_' + 'q1_' + str(q1Start) + '-' + str(q1Stop) + '_q2_' + str(q2Start) + '-' + str(q2Stop)+ '_'  + str(int(eta))
    unstableFileName = 'unstable_' + str(int(nIons)) + '_ions_' + str(int(nElectrons)) + '_electrons_' + 'q1_' + str(q1Start) + '-' + str(q1Stop) + '_q2_' + str(q2Start) + '-' + str(q2Stop)+ '_' + str(int(eta))
        
    unstableTriangles = []
    
    try:
        with open(r"data/triangles/" + unstableFileName + ".dat") as csvFile:
            csvReader = csv.reader(csvFile, delimiter = "\t")
                    
            for row in csvReader:
                
                triangle = []
                data = row[0].split(',')
                for i in range(3):
                    x = float(data[2*i])
                    y = float(data[2*i + 1])
                    triangle.append(np.array([x,y]))
                triangle = np.array([triangle, int(data[-1])],dtype=object)
                unstableTriangles.append(triangle)
            unstableTriangles = np.array(unstableTriangles)
    except Exception:
        pass
        
    stableTriangles = []
    
    try:    
        with open(r"data/triangles/" + stableFileName + ".dat") as csvFile:
            csvReader = csv.reader(csvFile, delimiter = "\t")
                    
            for row in csvReader:
                
                triangle = []
                data = row[0].split(',')
                for i in range(3):
                    x = float(data[2*i])
                    y = float(data[2*i + 1])
                    triangle.append(np.array([x,y]))
                triangle = np.array([triangle, int(data[-1])],dtype=object)
                stableTriangles.append(triangle)
            stableTriangles = np.array(unstableTriangles)  
    except Exception:
        pass

    return unstableTriangles, stableTriangles