import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
from matplotlib import cm           # Colour maps for the contour graph

from integration import *
from createSystem import *
from timeit import default_timer as timer
from multiprocessing import Pool

def IsStable(stabilityParam):
    """    if stabilityParam < 0.8:
            return 0
        else:
            return 1
    """
    return stabilityParam

def TriangleTest(p1, p2, p3):
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

def InTriangle(p, regions):
    
    for region in regions:
        triangle, inValue = region
        p1, p2, p3 = triangle
    
        inCheck1 = TriangleTest(p, p1, p2) < 0.0
        inCheck2 = TriangleTest(p, p2, p3) < 0.0
        inCheck3 = TriangleTest(p, p3, p1) < 0.0
        
        if((inCheck1 == inCheck2) and (inCheck2 == inCheck3)):
            return True
    
    return False    

def MakePoolList(system, ODESystem, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, tmax, dt):
    
    q1Step = (q1Stop - q1Start) / q1Resol
    q2Step = (q2Stop - q2Start) / q2Resol
    
    n = len(system)
    m = 0    
    for i in range(n):
        if(system[i,3] > 0):
            m = m + 1
    nParticles = (m, n-m)#numer of (ions, electrons)

    loadParams = np.array([q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, int(f2/f1), None, None],dtype=object)#need to get all this data    
    unstableRegion, stableRegion = LoadTriangles(loadParams)
    
    args = []
    
    a = 0    
    q1 = q1Start
    for i in range(q1Resol):
        q2 = q2Start
        for j in range(q2Resol):
            
            p = np.array([q2, q1])
            if InTriangle(p, unstableRegion):
                stabilityValue = n
            elif InTriangle(p, stableRegion):
                stabilityValue = 0
            else:
                stabilityValue = -1
                
            params = [system, np.array([a,q1,q2]), tmax, dt, ODESystem, loadParams]#loadParams are needed for ditching regions(parameter for LoadTriangles)
            args.append(tuple([params, i, j, stabilityValue]))
            
            q2 = q2 + q2Step
        q1 = q1 + q1Step
        
    return np.array(args, dtype=object)

def IntWrapper(params, i, j, stabilityValue):
    
    allowDitching = True
    
    system, trapParams, tmax, dt, ODESystem, loadParams = params
    n = len(system)
    
    if allowDitching:
        
        if stabilityValue == -1:
            stabilityValue = ODEint(system, trapParams, tmax, dt, ODESystem)[-1]
            
    else:
        stabilityValue = ODEint(system, trapParams, tmax, dt, ODESystem)[-1]
        
    result = np.array([IsStable(stabilityValue), i, j])
    return result

def StabilityDiagram(system, ODESystem, q1Start=0.0, q1Stop=0.15, q1Resol=20, q2Start=0.0, q2Stop=1.0, q2Resol=20, tmax=1.3e+2, dt=1e-2):
    
    stability = np.zeros((q1Resol,q2Resol))

    start = timer()

    pool = Pool()# take maximum available number of cpus
    args = MakePoolList(system, ODESystem, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, tmax, dt)
    results = pool.starmap(IntWrapper, args)
    
    for result in results:
        stabilityValue, i, j = result
        stability[i,j] = stabilityValue
        
    pool.close()
    stop = timer()
    time = stop-start
    print(time)
    
    n = len(system)
    m = 0    
    for i in range(n):
        if(system[i,3] > 0):
            m = m + 1
    nParticles = (m, n-m)#numer of (ions, electrons)
        
    params = np.array([q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, time, f1, f2], dtype=object)
    SaveStabilityDiagram(stability, params)
  
    return stability, params   

triangleStable = []
triangleUnstable = []
def ClickStabilityRegions(fileName='0_ions_1_electrons_q1_0-0.06_q2_0-0.48_700x700_13'):    
    data, params = LoadStabilityDiagram(fileName)
    
    def SaveTriang(tStable=triangleStable, tUnstable=triangleUnstable, parameters=params):
        SaveTriangles(tUnstable, tStable, parameters)  
        
    q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, eta, f1, f2 = params
    
    fig = plt.figure()
        
    x_vals = np.linspace(q2Start, q2Stop, int(q2Resol))
    y_vals = np.linspace(q1Start, q1Stop, int(q1Resol))
    X, Y = np.meshgrid(x_vals, y_vals)  

    plt.xlabel('$q_{2}$')
    plt.ylabel('$q_{1}$')
    plt.contourf(X, Y, data)
                
    def onclick(event):
        global triangleStable
        global triangleUnstable
        
        if not((event.button is MouseButton.RIGHT)or(event.button is MouseButton.LEFT)):
            fig.canvas.mpl_disconnect(cid)
            
        x, y = event.xdata, event.ydata
        if(x != None):
            print(x,y)
            p = [x,y]
                
            if event.button is MouseButton.LEFT:
                triangleUnstable.append(p)
                
            if event.button is MouseButton.RIGHT:
                triangleStable.append(p)

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    
    return params
##############################testing
def MakePoolListEdge(system, ODESystem, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, tmax, dt, needComputationPack):
    
    q1Step = (q1Stop - q1Start) / q1Resol
    q2Step = (q2Stop - q2Start) / q2Resol
    
    needComputation, previous_q1Resol, previous_q2Resol = needComputationPack
    
    n = len(system)
    m = 0    
    for i in range(n):
        if(system[i,3] > 0):
            m = m + 1
    nParticles = (m, n-m)#numer of (ions, electrons)

    loadParams = np.array([q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, int(f2/f1), None, None],dtype=object)#need to get all this data    
    
    args = []
    
    a = 0    
    q1 = q1Start
    for i in range(q1Resol):
        q2 = q2Start
        for j in range(q2Resol):
            previous_i = int(np.round(i * previous_q1Resol / q1Resol))            
            previous_j = int(np.round(j * previous_q2Resol / q2Resol))            
            
            stabilityValue = needComputation[previous_i, previous_j]
                
            params = [system, np.array([a,q1,q2]), tmax, dt, ODESystem, loadParams]#loadParams are needed for ditching regions(parameter for LoadTriangles)
            args.append(tuple([params, i, j, stabilityValue]))
            
            q2 = q2 + q2Step
        q1 = q1 + q1Step
        
    return np.array(args, dtype=object)

def IntWrapperEdge(params, i, j, stabilityValue):
    
    allowDitching = True
    
    system, trapParams, tmax, dt, ODESystem, loadParams = params
    n = len(system)
    
    if allowDitching:
        
        if stabilityValue == -1:
            stabilityValue = ODEint(system, trapParams, tmax, dt, ODESystem)[-1]
            
    else:
        stabilityValue = ODEint(system, trapParams, tmax, dt, ODESystem)[-1]
        
    result = np.array([IsStable(stabilityValue), i, j])
    return result

def StabilityDiagramEdge(system, ODESystem, previousFile, q1Start=0.0, q1Stop=0.15, q1Resol=20, q2Start=0.0, q2Stop=1.0, q2Resol=20, tmax=1.3e+2, dt=1e-2):
    
    stability = np.zeros((q1Resol, q2Resol))
    needComputationPack = PrepForStabilityEdge(previousFile)

    start = timer()

    pool = Pool()# take maximum available number of cpus
    args = MakePoolListEdge(system, ODESystem, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, tmax, dt, needComputationPack)
    results = pool.starmap(IntWrapper, args)
    
    for result in results:
        stabilityValue, i, j = result
        i, j = int(i), int(j)
        stability[i,j] = stabilityValue
        
    pool.close()
    stop = timer()
    time = stop-start
    print(time)
    
    n = len(system)
    m = 0    
    for i in range(n):
        if(system[i,3] > 0):
            m = m + 1
    nParticles = (m, n-m)#numer of (ions, electrons)
        
    params = np.array([q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, time, f1, f2], dtype=object)
    SaveStabilityDiagram(stability, params)
  
    return stability, params   


def PrepForStabilityEdge(fileName='0_ions_1_electrons_q1_0-0.06_q2_0-0.48_700x700_13'):
    
    stability, params = LoadStabilityDiagram(fileName)
    q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, eta, f1, f2 = params
    nIons, nElectrons = nParticles
    n = nIons + nElectrons
    
    q1Step = (q1Stop - q1Start) / q1Resol
    q2Step = (q2Stop - q2Start) / q2Resol
    
    needComputation = copy(stability)
    
    def TalkToRightNeighbour(i, j):
        if stability[i,j] != stability[i,j+1]:
            needComputation[i,j] = -1
            needComputation[i,j+1] = -1
            
    def TalkToTopNeighbour(i, j):
        if stability[i,j] != stability[i+1,j]:
            needComputation[i,j] = -1
            needComputation[i+1,j] = -1
        
    def GoToNeighbour(i,j):
        if not((i == q1Resol) and (j == q2Resol)):
            if i == q1Resol:
                TalkToRightNeighbour(i,j)
                
            elif j == q2Resol:
                TalkToTopNeighbour(i,j)
                
            else:
                TalkToRightNeighbour(i,j)
                TalkToTopNeighbour(i,j) 
               
    for i in range(q1Resol-1):
        for j in range(q2Resol-1):
            GoToNeighbour(i, j)
            
    return np.array([needComputation, q1Resol, q2Resol], dtype=object)
