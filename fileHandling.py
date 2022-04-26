"""
In this file we handle file input and output
"""
import numpy as np
from parameters import *
from plotting import *
from integration import *
from copy import copy
import csv

"""
I've decided to put some information about the system into the corresponding file name. 
Function ParseFileName() extracts this information. Its single argument is the name of the file
we want to get information about.
"""
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


"""
Function MakeFileName() returns a string, which is a name of a file in a format compatible with the function ParseFileName().
We set its second argument withResolution to False if we don't want to have resolution in the name of the file.
We use this option when computing stability diagram just on the edge of the stability region.
"""
def MakeFileName(params, withResolution=True):
    
    q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, time, f1, f2 = params
    nIons, nElectrons = nParticles
    
    if((f1 == None)or(f2 == None)):
        f2 = 1337
        f1 = 1
    
    if withResolution:
        fileName = str(int(nIons)) + '_ions_' + str(int(nElectrons)) + '_electrons_' + 'q1_' + str(q1Start) + '-' + str(q1Stop) + '_q2_' + str(q2Start) + '-' + str(q2Stop) + '_' + str(int(q2Resol)) + 'x' + str(int(q1Resol)) + '_' + str(int(f2/f1))
    else:
        fileName = str(int(nIons)) + '_ions_' + str(int(nElectrons)) + '_electrons_' + 'q1_' + str(q1Start) + '-' + str(q1Stop) + '_q2_' + str(q2Start) + '-' + str(q2Stop) + '_' + str(int(f2/f1))

    
    return fileName

"""
SaveParticleSystem() is a function that saves position, velocity, mass and charge of all particles in the given system into a file.
It has two arguments: the system we want to save and the name of the file we want it to save to.
"""
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
            
"""
LoadParticleSystem() is a function which returns particle system from a file. 
The file must be prepared in the format: position, velocity, mass, charge -> same as file created by function SaveParticleSystem().
"""
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

"""
SaveStabilityDiagram() saves stability matrix into a file
"""            
def SaveStabilityDiagram(stability, params):
    
    fileName = MakeFileName(params)
    
    with open(r"data/stability_diagrams/" + fileName + ".dat","w", newline="") as csvFile:
        csvWriter = csv.writer(csvFile, delimiter = "\t")
        #csvWriter.writerow("Ta toto su nejake data")
        for row in stability:
            csvWriter.writerow([str(row)[1:-1]])
            
"""
LoadStabilityDiagram() loads stability matrix from a file
"""            
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
            
    params = np.array([q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, np.array([nIons, nElectrons]), eta, f1, f2], dtype=object)

    return stability, params
"""
To reduce time needed for making a diagram of stability, we can choose the regions in which we are certain stability (or instability)
of the solution. Than we save those regions (in shape of triagles) by the function SaveTriangles() into the file.
It has three arguments: first -> triangleUnstable is the list all triangles in unstable region along with the actual value of stability
(value of stability is needed for stability diagram of the system with multiple particles, as it represents number of unstable particles)
"""    
def SaveTriangles(triangleUnstable, triangleStable, params):
    
    q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, eta, f1, f2 = params
    nIons, nElectrons = nParticles
    n = nIons + nElectrons
    
    nUnstable = len(triangleUnstable) // 3
    unstables = []
    for i in range(nUnstable):
        unstables.append([triangleUnstable[i*3:i*3 + 3], n])#saves 3 points unambiguously defining a triangle and saving number of unstable particles aswell
        
    
    nStable = len(triangleStable) // 3
    stables = []
    for i in range(nStable):
        stables.append([triangleStable[i*3:i*3 + 3], 0])#saves 3 points unambiguously defining a triangle and saving number of unstable particles which is zero
    
    fileNameStable = 'stable_' + MakeFileName(params, withResolution=False)
    fileNameUnstable = 'unstable_' + MakeFileName(params, withResolution=False)
    
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
            
"""
Function LoadTriangles() returns unstable and stable regions loaded from a file
"""            
def LoadTriangles(params):

    fileNameStable = 'stable_' + MakeFileName(params, withResolution=False)
    fileNameUnstable = 'unstable_' + MakeFileName(params, withResolution=False)
        
    unstableTriangles = []
    
    try:
        with open(r"data/triangles/" + fileNameUnstable + ".dat") as csvFile:
            csvReader = csv.reader(csvFile, delimiter = "\t")
                    
            for row in csvReader:
                
                triangle = []
                data = row[0].split(',')
                for i in range(3):
                    x = float(data[2*i])
                    y = float(data[2*i + 1])
                    triangle.append(np.array([x,y]))
                triangle = np.array([triangle, int(data[6])],dtype=object)
                unstableTriangles.append(triangle)
            unstableTriangles = np.array(unstableTriangles)
    except Exception:
        pass
        
    stableTriangles = []
    
    try:    
        with open(r"data/triangles/" + fileNameStable + ".dat") as csvFile:
            csvReader = csv.reader(csvFile, delimiter = "\t")
                    
            for row in csvReader:
                
                triangle = []
                data = row[0].split(',')
                for i in range(3):
                    x = float(data[2*i])
                    y = float(data[2*i + 1])
                    triangle.append(np.array([x,y]))
                triangle = np.array([triangle, int(data[6])],dtype=object)
                stableTriangles.append(triangle)
            stableTriangles = np.array(stableTriangles)  
    except Exception:
        pass

    return unstableTriangles, stableTriangles