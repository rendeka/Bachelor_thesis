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

def InTriangle(p, triangles):
    
    for triangle in triangles:
        p1, p2, p3, inValue = triangle
    
        inCheck1 = TriangleTest(p, p1, p2) < 0.0
        inCheck2 = TriangleTest(p, p2, p3) < 0.0
        inCheck3 = TriangleTest(p, p3, p1) < 0.0
        
        if((inCheck1 == inCheck2) and (inCheck2 == inCheck3)):
            return True
    
    return False

"""
def MakeTriangles(data, params):
    
    q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, time, f1, f2 = params
    
    fig = plt.figure()
        
    x_vals = np.linspace(q2Start, q2Stop, int(q2Resol))
    y_vals = np.linspace(q1Start, q1Stop, int(q1Resol))
    X, Y = np.meshgrid(x_vals, y_vals) 

    plt.xlabel('$q_{2}$')
    plt.ylabel('$q_{1}$')
    plt.contourf(X, Y, data)
        
    def onclick(event):
        if not((event.button is MouseButton.RIGHT)or(event.button is MouseButton.LEFT)):
            fig.canvas.mpl_disconnect(cid)
        x, y = event.xdata, event.ydata
        print(x,y)
        p = np.array([x,y])
        
        if event.button is MouseButton.RIGHT:
            print('unstable')
            
        if event.button is MouseButton.LEFT:
            print('stable')
        global coordsC
        coordsC.append(p)

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
"""
        
def ComputationNeeded(trapParams, n):#we dont want to waste time computing stability in regions far from the edge of stability

    a, q1, q2 = trapParams
    p = np.array([q2, q1])

    triangles = []
    
    #""" #for f2/f1 = 13
    p1 = np.array([0,0.02])
    p2 = np.array([0,0.14])
    p3 = np.array([0.6,0.14])    
    triangles.append(np.array([p1, p2, p3, n],dtype=object))# n for instable region
    
    p1 = np.array([0.8,0.14])
    p2 = np.array([0.9,0.10])
    p3 = np.array([0.9,0.14])    
    triangles.append(np.array([p1, p2, p3, n],dtype=object))
    
    p1 = np.array([0.3,0])
    p2 = np.array([0.8,0])
    p3 = np.array([0.65,0.06])    
    triangles.append(np.array([p1, p2, p3, 0],dtype=object))# 0 for stable region
    #"""
    
    """# for f2/f1 = 170
    p1 = np.array([0.47,0.04])
    p2 = np.array([0.41,0.04])
    p3 = np.array([0.47,0])    
    triangles.append(np.array([p1, p2, p3, n],dtype=object))# n for instable region

    p1 = np.array([0,0])
    p2 = np.array([0.4,0.04])
    p3 = np.array([0,0.04])    
    triangles.append(np.array([p1, p2, p3, n],dtype=object))# n for instable region
  
    p1 = np.array([0.15,0])
    p2 = np.array([0.36,0.025])
    p3 = np.array([0.43,0])    
    triangles.append(np.array([p1, p2, p3, 0],dtype=object))
    """    
    
    """# for f2/f1 = 13 but asymetric equation
    p1 = np.array([0,0.01])
    p2 = np.array([0.4,0.06])
    p3 = np.array([0,0.06])    
    triangles.append(np.array([p1, p2, p3, n],dtype=object))# n for instable region

    p1 = np.array([0.15,0])
    p2 = np.array([0.4,0.03])
    p3 = np.array([0.42,0])    
    triangles.append(np.array([p1, p2, p3, 0],dtype=object))
    """    
    
    """# for f2/f1 = 13 but asymetric equation
    p1 = np.array([0,0.003])
    p2 = np.array([0,0.045])
    p3 = np.array([0.4,0.045])    
    triangles.append(np.array([p1, p2, p3, n],dtype=object))# n for instable region

    p1 = np.array([0.2,0])
    p2 = np.array([0.5,0])
    p3 = np.array([0.38,0.025])    
    triangles.append(np.array([p1, p2, p3, 0],dtype=object))
    #"""  
    
    triangles = np.array(triangles)
    
    return InTriangle(p, triangles)
    
"""
def IntWrapper(params, i, j):
    
    allowDitching = True
    
    system, trapParams, tmax, dt, ODESystem = params
    
    if allowDitching:
        subResult = ComputationNeeded(trapParams, len(system))
        if(subResult == -1):
            stabilityValue = ODEint(system, trapParams, tmax, dt, ODESystem)[-1]
        else:
            stabilityValue = subResult
    else:
        stabilityValue = ODEint(system, trapParams, tmax, dt, ODESystem)[-1]
        
    result = np.array([IsStable(stabilityValue), i, j])
    return result
"""

def MakePoolList(system, ODESystem, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, tmax, dt):
    
    q1Step = (q1Stop - q1Start) / q1Resol
    q2Step = (q2Stop - q2Start) / q2Resol
    
    n = len(system)
    m = 0
    
    for i in range(n):
        if(system[i,3] > 0):
            m = m + 1
    nParticles = (m, n-m)#numer of ions and electrons
    
    a = 0
    
    args = []
    loadParams = np.array([q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, int(f2/f1), None, None],dtype=object)#need to get all this data
    
    q1 = q1Start
    for i in range(q1Resol):
        q2 = q2Start
        for j in range(q2Resol):
            params = [system, np.array([a,q1,q2]), tmax, dt, ODESystem, loadParams]#loadParams are needed for ditching regions(parameter for LoadTriangles)
            args.append(tuple([params, i, j]))
            q2 = q2 + q2Step
        q1 = q1 + q1Step
        
    return np.array(args, dtype=object)

def IntWrapper(params, i, j):
    
    allowDitching = True
    
    system, trapParams, tmax, dt, ODESystem, loadParams = params
    n = len(system)
    
    if allowDitching:
        
        a, q1, q2 = trapParams
        p = np.array([q2, q1])
                
        unstableRegion, stableRegion = LoadTriangles(loadParams)
        
        if InTriangle(p, unstableRegion):
            stabilityValue = n
        elif InTriangle(p, stableRegion):
            stabilityValue = 0
        else:
            stabilityValue = ODEint(system, trapParams, tmax, dt, ODESystem)[-1]

    else:
        stabilityValue = ODEint(system, trapParams, tmax, dt, ODESystem)[-1]
        
    result = np.array([IsStable(stabilityValue), i, j])
    return result

def StabilityDiagram(system, ODESystem, q1Start=0, q1Stop=0.15, q1Resol=20, q2Start=0, q2Stop=1, q2Resol=20, tmax=1.3e+2, dt=1e-2):
    
    stability = np.zeros((q2Resol,q1Resol))

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
    nParticles = (m, n-m)#numer of ions and electrons
        
    params = np.array([q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, time, f1, f2], dtype=object)
    SaveStabilityDiagram(stability, params)
  
    return stability, params   
