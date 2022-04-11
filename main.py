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

def MakeCoulombCrystal():
    system, solution = IonCrystal(1)
    
    while len(system)<26:
        
        #input("Press Enter to continue...")
        print(len(system))
        system = AddIon(system)
        system, solution = IonCrystal(system=system)
        
    system, solution = IonCrystal(system=system)
    
    SaveParticleSystem(system, 'coulomb_crystals/crystal2')
        
    PlotODESolution(solution)
    PlotFinalPositions(solution)
    
    
def TestCoulombCrystal():
    system = LoadParticleSystem('coulomb_crystals/crystal')
    trapParams = np.array([0, 0.4 * ionMass / electronMass, 0])
    
    for particle in system:
        particle[1] = np.zeros(3)
        
    solution = ODEint(system, trapParams, endTime, timeStep, ODESystem=ODESystemExact,  Step=StepEulerAdvanced)      
    #PlotODESolution(solution)
    #PlotFinalPositions(solution)
    
    SaveParticleSystem(system, 'coulomb_crystals/crystal-evolution')
    
def PlotCrystalTest():
    s=LoadParticleSystem('coulomb_crystals/crystal-evolution') 
    sol = ODEint(s, np.array([0,0,0]), timeStep, timeStep, ODESystem=ODESystemExact,  Step=StepEulerAdvanced)
    PlotFinalPositions(sol)
    
def SolveParticleSystem():
    initsystem = MakeParticleSystem(0,1)
    SaveParticleSystem(initsystem, '1')

    #initsystem = LoadParticleSystem('1') 
    #initsystem = LoadParticleSystem('oneElectron')
    
    #trapParams = np.array([0, 0.12, 0.45])
    trapParams = np.array([0, 0.05, 0.45])#test
    
    ODESolution = ODEint(initsystem, trapParams, 2*endTime, timeStep, ODESystemExact, StepEulerAdvanced)
    
    #PlotCoordinates(ODESolution)
    PlotEnergy(ODESolution)
    PlotODESolution(ODESolution)

    
def MakeStabilityDiagram():
    initsystem = LoadParticleSystem('1')
    #initsystem = MakeParticleSystem(0,1)
    #SaveParticleSystem(initsystem, '2')
    stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, 0, 0.14, 10, 0, 0.9, 10, endTime, timeStep) 

    PlotStability(stability, plotParams)

if __name__ == '__main__':
    
    timeStep = 1/(200) #in normal time scale 1/(100 * f2) -> faster period divided into 100 segments 
    endTime = 5 * f2 / f1 #in normal time scale 5 * (1 / f1) -> five slower periods

    #MakeCoulombCrystal() 
    
    #TestCoulombCrystal()
    
    #PlotCrystalTest() 
    
    #SolveParticleSystem()
    
    #MakeStabilityDiagram()
