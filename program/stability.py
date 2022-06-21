"""
In this module we make diagrams of stability depending on trap parameters setting.
We also use two approaches to try to decrease time needed for making such a diagram.
First approach is to make low resolution picture first, then mark regions which you are sure to be stable
or unstable with help of function ClickStabilityRegions() and in this regions we won't solve equations of motion
for much higher resolution, saving precious time.
The second approach is to gradually compute stability diagram for higher and higher resolutions, but in this case we
solve equations of motion only on the edge between stable and unstable regions, which increase the speed of computation radically.
(Functions working this way have the name ending with 'Edge'). This approach is much faster but has it's own drawbacks.
If we start with the picture with low resolution and then continue to improove quality just od the edge of stability,
then we won't by able to detect thin stability strips that escaped our sight in low resolution setting.
"""

"""
stability for light particle is stable for heavy particle if the condition q1 * (f2 / f1)**2 * (m2 / m1) < 0.9
"""

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
from matplotlib import cm           # Colour maps for the contour graph

from integration import *
#from createSystem import *
from fileHandling import *
#from parameters import *

from scipy.sparse.linalg import eigsh


from timeit import default_timer as timer
from multiprocessing import Pool

#from scipy.sparse.linalg import eigsh


def IsStable(stabilityParam):
    """
    Function IsStable() decides whether given stability parameters corresponds to stable or unstable trajectory.
    In this case it doesn't do anything.
    """
    result = stabilityParam
    return result


def TriangleTest(p1, p2, p3):
    """
    Function TriangleTest() is given 3 points in x-y plane. Points p2 and p3 define a straight line dividing x-y plane into two.
    TriangleTest() checks on which half-plane the point p1 is. 
    """
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])


def InTriangle(p, regions):
    """
    Function InTriangle() has two arguments. First is p -> position in x-y plane. Second is list of tringles.
    If point p belongs in some of the given triangles, InTriangle() returns True, else it returns False
    """
    
    for region in regions:
        triangle, inValue = region
        p1, p2, p3 = triangle
    
        inCheck1 = TriangleTest(p, p1, p2) < 0.0
        inCheck2 = TriangleTest(p, p2, p3) < 0.0
        inCheck3 = TriangleTest(p, p3, p1) < 0.0
        
        if((inCheck1 == inCheck2) and (inCheck2 == inCheck3)):
            return True
    
    return False    

def MakePoolList(system, ODESystem, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, tmax, dt, freezeIons, velocityDiagram):
    """
    Prepares list of arguments for parallel computing. Each element from this list is
    argument for function IntWrapper().
    """
    
    q1Step = (q1Stop - q1Start) / q1Resol
    q2Step = (q2Stop - q2Start) / q2Resol
    
    n = len(system)
    m = 0    
    for i in range(n):
        if(system[i,3] > 0):
            m = m + 1
    nParticles = (m, n-m)#numer of (ions, electrons)

    loadParams = np.array([q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, int(f2/f1), f1, f2],dtype=object)  
    unstableRegion, stableRegion = LoadTriangles(loadParams)
    
    args = []
    
    a = 0    
    q1 = q1Start
    for i in range(q1Resol):
        q2 = q2Start
        for j in range(q2Resol):
            
            if velocityDiagram:
                stabilityValue = -1
            
            else:
                p = np.array([q2, q1])
                if InTriangle(p, unstableRegion):# here we decide whether we need to compute the stability value or whether we already know it for given parameters
                    stabilityValue = n
                elif InTriangle(p, stableRegion):
                    stabilityValue = 0
                else:
                    stabilityValue = -1
                
            params = [system, np.array([a,q1,q2]), tmax, dt, ODESystem, freezeIons, velocityDiagram]
            args.append(tuple([params, i, j, stabilityValue]))
            
            q2 = q2 + q2Step
        q1 = q1 + q1Step
        
    return np.array(args, dtype=object)

def IntWrapper(params, i, j, stabilityValue):
    """
    Function IntWrapper returns stability value for given set of trapping paramaters.
    """
        
    system, trapParams, tmax, dt, ODESystem, freezeIons, velocityDiagram = params
    n = len(system)
            
    if stabilityValue == -1: #stability value is computed only if 
        stabilityValue = ODEint(system, trapParams, tmax, dt, ODESystem, freezeIons=freezeIons, velocityDiagram=velocityDiagram)[-1]
            
        
    result = np.array([IsStable(stabilityValue), i, j])
    return result

def StabilityDiagram(system, ODESystem, q1Start=0.0, q1Stop=0.15, q1Resol=20, q2Start=0.0, q2Stop=1.0, q2Resol=20, tmax=1.3e+2, dt=1e-2, freezeIons=False, velocityDiagram=False):
    
    stability = np.zeros((q1Resol,q2Resol))

    start = timer()

    pool = Pool()# take maximum available number of cpus
    args = MakePoolList(system, ODESystem, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, tmax, dt, freezeIons, velocityDiagram)
    results = pool.starmap(IntWrapper, args)
    
    for result in results:
        stabilityValue, i, j = result
        i, j = int(i), int(j)
        stability[i,j] = stabilityValue
        
    pool.close()
    stop = timer()
    time = stop-start
    print('time: ', time, 'seconds')
    
    n = len(system)
    m = 0    
    for i in range(n):
        if(system[i,3] >= 0):
            m = m + 1
    nParticles = (m, n-m)#numer of (ions, electrons)
        
    params = np.array([q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, time, f1, f2], dtype=object)
    SaveStabilityDiagram(stability, params, velocityDiagram=velocityDiagram)
  
    return stability, params   

triangleStable = []
triangleUnstable = []
def ClickStabilityRegions(fileName='0_ions_1_electrons_q1_0-0.06_q2_0-0.48_700x700_13'):    
    data, params = LoadStabilityDiagram(fileName)
    
    def SaveTriang(tStable=triangleStable, tUnstable=triangleUnstable, parameters=params):
        SaveTriangles(tUnstable, tStable, parameters)  
        
    if fileName[0] == 'd':
        q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, eta = ParseFileNameDet(fileName)
    else:
        q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nIons, nElectrons, eta = ParseFileName(fileName)

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

"""Here follows similar function as we've already seen, but for computing stability just od the edge of stability"""

def MakePoolListEdge(system, ODESystem, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, tmax, dt, needComputationPack, freezeIons):
    
    q1Step = (q1Stop - q1Start) / q1Resol
    q2Step = (q2Stop - q2Start) / q2Resol
    
    needComputation, previous_q1Resol, previous_q2Resol = needComputationPack
    
    previous_q1Step = (q1Stop - q1Start) / previous_q1Resol
    previous_q2Step = (q2Stop - q2Start) / previous_q2Resol
    
    n = len(system)
    m = 0    
    for i in range(n):
        if(system[i,3] > 0):
            m = m + 1
    nParticles = (m, n-m)#number of (ions, electrons)
    
    args = []
    
    a = 0  
    
    previous_i = 0
    previous_j = 0
    
    q1 = q1Start
    for i in range(q1Resol):
        q2 = q2Start
        for j in range(q2Resol):
            previous_i = int(np.round(i * q1Step / previous_q1Step))            
            previous_j = int(np.round(j * q2Step / previous_q2Step))

            if previous_i == previous_q1Resol: previous_i = previous_q1Resol - 1 #it's obvious why we need to do this if you draw a picture..    
            if previous_j == previous_q2Resol: previous_j = previous_q2Resol - 1         
            
            stabilityValue = needComputation[previous_i, previous_j]
                
            params = [system, np.array([a,q1,q2]), tmax, dt, ODESystem, freezeIons]
            args.append(tuple([params, i, j, stabilityValue]))
            
            q2 = q2 + q2Step
        q1 = q1 + q1Step
        
    return np.array(args, dtype=object)

def IntWrapperEdge(params, i, j, stabilityValue):
        
    system, trapParams, tmax, dt, ODESystem, freezeIons = params
    n = len(system)
            
    if stabilityValue == -1:
        stabilityValue = ODEint(system, trapParams, tmax, dt, ODESystem, freezeIons=freezeIons)[-1]
        
    result = np.array([IsStable(stabilityValue), i, j])
    return result

def StabilityDiagramEdge(system, ODESystem, previousFile, q1Start=0.0, q1Stop=0.15, q1Resol=20, q2Start=0.0, q2Stop=1.0, q2Resol=20, tmax=1.3e+2, dt=1e-2, freezeIons=False):
    
    stability = np.zeros((q1Resol, q2Resol))
    needComputationPack = PrepForStabilityEdge(previousFile)
    
    pool = Pool()# take maximum available number of cpus
    args = MakePoolListEdge(system, ODESystem, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, tmax, dt, needComputationPack, freezeIons)
    results = pool.starmap(IntWrapperEdge, args)
    
    for result in results:
        stabilityValue, i, j = result
        i, j = int(i), int(j)
        stability[i,j] = stabilityValue
        
    pool.close()
    
    n = len(system)
    m = 0    
    for i in range(n):
        if(system[i,3] > 0):
            m = m + 1
    nParticles = (m, n-m)#numer of (ions, electrons)
        
    params = np.array([q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, None, f1, f2], dtype=object)
    SaveStabilityDiagram(stability, params)
  
    return stability, params   


def PrepForStabilityEdge(fileName, needComputationValue=-1):
    """
    Function PrepForStabilityEdge() loads the stability diagram which is represented by AxB matrix filled with stability values.
    It returns same AxB matrix but if some value of stability isn't the same as values of it's closest neighbours,
    then it changes these stability values to -1. Which means that we will compute the equations of motion in thiw regions
    with the better resolution in next interation.
    """
    
    stability, params = LoadStabilityDiagram(fileName)
    q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, eta, f1, f2 = params
    nIons, nElectrons = nParticles
    n = nIons + nElectrons
        
    q1Step = (q1Stop - q1Start) / q1Resol
    q2Step = (q2Stop - q2Start) / q2Resol
    
    needComputation = copy(stability)
    
    def TalkToRightNeighbour(i, j):
        if stability[i,j] != stability[i,j+1]:
            needComputation[i,j] = needComputationValue
            needComputation[i,j+1] = needComputationValue
            
    def TalkToTopNeighbour(i, j):
        if stability[i,j] != stability[i+1,j]:
            needComputation[i,j] = needComputationValue
            needComputation[i+1,j] = needComputationValue
            
    def TalkToUpperDiagonalNeighbour(i, j):
        if stability[i,j] != stability[i+1,j+1]:
            needComputation[i,j] = needComputationValue
            needComputation[i+1,j+1] = needComputationValue
    
    def TalkToLowerDiagonalNeighbour(i, j):
        if stability[i,j] != stability[i-1,j+1]:
            needComputation[i,j] = needComputationValue
            needComputation[i-1,j+1] = needComputationValue
            
    
    
    """
    We need to check every pair of closest neighbours in a grid. To do so, we start in the bottom left corner and we check
    the neighbour to the right and then top neighbour. Then we repeat this process through all rows and columns(if there is no more neighbour
                                                                                                                to the right or top then we don't check it)
    The function GoToNeighbour() does exactly that
    """
    def GoToNeighbour(i,j):
        if not((i == q1Resol-1) and (j == q2Resol-1)):
            
            if i == q1Resol-1:
                TalkToRightNeighbour(i,j)
                TalkToLowerDiagonalNeighbour(i,j)                            
                
            elif j == q2Resol-1:
                TalkToTopNeighbour(i,j)  
            
            elif i == 0:
                TalkToRightNeighbour(i,j)
                TalkToTopNeighbour(i,j)
                TalkToUpperDiagonalNeighbour(i,j)
                
            else:
                TalkToRightNeighbour(i,j)
                TalkToTopNeighbour(i,j)
                TalkToUpperDiagonalNeighbour(i,j)
                TalkToLowerDiagonalNeighbour(i,j)
    
    """In next for loops we are making new AxB matrix with stability value = -1 on the edge of two regions with different stability"""
    for i in range(q1Resol):
        for j in range(q2Resol):
            GoToNeighbour(i,j)
            
    return np.array([needComputation, q1Resol, q2Resol], dtype=object)

"""
next we will create stability diagram using floquet theory
"""

def StabilityMatrix(trapParams, f1, f2):
    
    a, q1, q2 = trapParams
    
    gcd = np.gcd(int(f1),int(f2))
    m = int(f2 // gcd)
    n = int(f1 // gcd)
    
    eta = f2 / f1
    
    resol = 10*m + 1
    
    stabilityMatrix = np.zeros((resol,resol))
    stabilityFactor = -2 #1    
    
    for i in range(resol):
        k = i - (resol-1)/2
        stabilityMatrix[i,i] = a * eta**2 - ((k)/m)**2
        if i + 2*n < resol:            
            stabilityMatrix[i,i + 2*n] = -q1 * stabilityFactor
            stabilityMatrix[i + 2*n,i] = -q1 * stabilityFactor
            
            if i + 2*m < resol:        
                stabilityMatrix[i,i + 2*m] = -q2 * stabilityFactor
                stabilityMatrix[i + 2*m,i] = -q2 * stabilityFactor
        
    return stabilityMatrix
    

def MakePoolListDet(q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, f1, f2):
    """
    Prepares list of arguments for parallel computing. Each element from this list is
    argument for function IntWrapper().
    """
    
    q1Step = (q1Stop - q1Start) / q1Resol
    q2Step = (q2Stop - q2Start) / q2Resol
    
    eta = int(f2/f1)
    
    loadParams = q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, eta, f1, f2
    unstableRegion, stableRegion = LoadTriangles(loadParams)

    args = []
    
    a = 0    
    q1 = q1Start
    for i in range(q1Resol):
        q2 = q2Start
        for j in range(q2Resol):
            
            p = np.array([q2, q1])            
            if InTriangle(p, unstableRegion):# here we decide whether we need to compute the stability value or whether we already know it for given parameters
                stabilityValue = 1
            elif InTriangle(p, stableRegion):
                stabilityValue = -1
            else:
                stabilityValue = 0            
                
            trapParams = np.array([a,q1,q2])
            args.append(tuple([trapParams, i, j, stabilityValue, f1, f2]))
            
            q2 = q2 + q2Step
        q1 = q1 + q1Step
        
    return np.array(args, dtype=object)

def IntWrapperDet(trapParams, i, j, stabilityValue, f1, f2):
    """
    Function IntWrapper returns stability value for given set of trapping paramaters.
    """
    if stabilityValue == 0:    
        matrix = StabilityMatrix(trapParams, f1, f2)            
        stabilityValue = -np.linalg.slogdet(matrix)[0]           
        
    result = np.array([IsStable(stabilityValue), i, j])
    return result

def StabilityDet(q1Start=0.0, q1Stop=0.15, q1Resol=20, q2Start=0.0, q2Stop=1.0, q2Resol=20, f1=3e6*2*np.pi, f2=17*3e6*2*np.pi):
    
    stability = np.zeros((q1Resol,q2Resol))

    start = timer()

    pool = Pool()# take maximum available number of cpus
    args = MakePoolListDet(q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, f1, f2)
    results = pool.starmap(IntWrapperDet, args)
    
    for result in results:
        stabilityValue, i, j = result
        i, j = int(i), int(j)
        stability[i,j] = stabilityValue
        
    pool.close()
    stop = timer()
    time = stop-start
    print('time: ', time, 'seconds')
    
    params = [q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, int(f2/f1)]
    
    SaveStabilityDiagramDet(stability, params)
    from plotting import PlotStabilityDet
    PlotStabilityDet(stability, params)
  
    #return stability, params   