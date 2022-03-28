import numpy as np

from integration import *
from createSystem import *
from timeit import default_timer as timer
from multiprocessing import Pool

def MakeList(system, ODESystem, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, tmax, dt):
    
    q1Step = (q1Stop - q1Start) / q1Resol
    q2Step = (q2Stop - q2Start) / q2Resol
    
    a = 0
    
    args = []
    
    q1 = q1Start
    for i in range(q1Resol):
        q2 = q2Start
        for j in range(q2Resol):
            params = [system, np.array([a,q1,q2]), tmax, dt, ODESystem]
            args.append(tuple([params, i, j]))
            q2 = q2 + q2Step
        q1 = q1 + q1Step
        
    return np.array(args, dtype=object)

def IntWrapper(params, i, j):

    system, trapParams, tmax, dt, ODESystem = params
    stabilityValue = ODEint(system, trapParams, tmax, dt, ODESystem)[-1]
    result = np.array([stabilityValue, i, j])
    return result

def StabilityDiagram(system, ODESystem, q1Start=0, q1Stop=0.15, q1Resol=20, q2Start=0, q2Stop=1, q2Resol=20, tmax=1.3e+2, dt=1e-2):
    
    stability = np.zeros((q1Resol,q2Resol))

    start = timer()

    pool = Pool()
    args = MakeList(system, ODESystem, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, tmax, dt)
    results = pool.starmap(IntWrapper, args)
    
    for result in results:
        stabilityValue, i, j = result
        stability[i,j] = stabilityValue
        
    pool.close()
    stop = timer()
    print(stop-start)
    
    nParticles = len(system)
    
    params = np.array([q1Resol, q2Resol, nParticles])
    SaveStabilityDiagram(stability, params)
  
    return stability, np.array([q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles])    