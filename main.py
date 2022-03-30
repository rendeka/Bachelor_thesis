#import numpy as np
#import matplotlib.pyplot as plt
#from timeit import default_timer as timer
#from copy import copy
#import sys

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
    
    tmax = 2.5 * f2 / f1 #in normal time scale 5 * (1 / f1) -> five slower periods
    dt = 1/(200) #in normal time scale 1/(100 * f2) -> faster period divided into 100 segments 

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
        
    solution = ODEint(system, trapParams, tmax, dt, ODESystem=ODESystemCrystal,  Step=StepEulerAdvanced)      
    PlotODESolution(solution)
    PlotFinalPositions(solution)
    """   
    
    
    """
    initsystem = MakeParticleSystem(0,1)
    #initsystem = LoadParticleSystem('1') 
    
    trapParams = np.array([0, 0.024, 0.37])
    
    ODESolution = ODEint(initsystem, trapParams, tmax, dt, ODESystemExact, StepEulerAdvanced)
    
    #PlotCoordinates(ODESolution)
    PlotEnergy(ODESolution)
    PlotODESolution(ODESolution)
    """
    
    
    #"""
    #initsystem = LoadParticleSystem('1')
    initsystem = MakeParticleSystem(0,1)

    stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, 0, 0.14, 160, 0, 0.55, 300, tmax, dt) 

    PlotStability(stability, plotParams)
    #sys.modules[__name__].__dict__.clear()
    #"""