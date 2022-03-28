#import numpy as np
#import matplotlib.pyplot as plt
#from timeit import default_timer as timer
#from copy import copy
import sys

#from scipy import interpolate
#from multiprocessing import Pool, Process, Value

#from coulomb import * 
#from equations import * 
#from intMethods import * 
from plotting import * 
#from parameters import * 
from integration import *
from createSystem import *
from stability import *


if __name__ == '__main__':
    """
    system = IonCrystal(1)
    
    while len(system)<26:
        
        #input("Press Enter to continue...")
        print(len(system))
        system = AddIon(system)
        system = IonCrystal(system=system)
    """ 
    
    
    """    
    system = LoadParticleSystem('coulomb_crystal_16')
    trapParams = [0, 4.5, 0]
    
    for particle in system:
        particle[1] = np.zeros(3)
        
    solution = ODEint(system, trapParams, tmax=1e-1, dt=1e-3, ODESystem=ODESystemCrystal,  Step=StepEulerAdvanced)      
    PlotODESolution(solution)
    PlotFinalPositions(solution)
    """   
    
    
    """
    initsystem = MakeParticleSystem(5,0)
    #initsystem = ReadParticleSystem()
    
    timeStep = 1e-2
    endTime = 5e+2    
    
    trapParams = np.array([0, 0, 0.5])
    
    ODESolution = ODEint(initsystem, trapParams, endTime, timeStep, ODESystemExact, StepEulerAdvanced)
    
    #PlotCoordinates(ODESolution)
    #PlotEnergy(ODESolution)
    PlotODESolution(ODESolution)
    """
    
    
    #"""
    initsystem = MakeParticleSystem(6,0)
    stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, 0, 0.01, 30, 0.5, 1, 30, 5e+1, 1e-2) 
    PlotStability(stability, plotParams)
    sys.modules[__name__].__dict__.clear()
    #"""