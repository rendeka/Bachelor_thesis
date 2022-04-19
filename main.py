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

def MakeCoulombCrystal(nCrystal = 36):
    system, solution = IonCrystal(1)
    
    while len(system) < nCrystal:
        
        #input("Press Enter to continue...")
        print(len(system))
        system = AddIon(system)
        system, solution = IonCrystal(system=system)
        
    system, solution = IonCrystal(system=system)
    
    SaveParticleSystem(system, 'coulomb_crystals/' + str(int(nCrystal)))
        
    PlotODESolution(solution)
    PlotFinalPositions(solution)
    
    
def TestCoulombCrystal(nCrystal='36'):
    system = LoadParticleSystem('coulomb_crystals/' + nCrystal)
    trapParams = np.array([0, 0.4, 0])
    
    for particle in system:
        particle[1] = np.zeros(3)
        
    solution = ODEint(system, trapParams, 1*endTime, timeStep, ODESystem=ODESystemExact,  Step=StepEulerAdvanced)      
    PlotODESolution(solution)
    PlotFinalPositions(solution)
    print(TotalVelocity(system))
    
    SaveParticleSystem(system, 'coulomb_crystals/crystal-evolution_' + nCrystal)
    
def PlotCrystalTest(nCrystal='36'):
    s = LoadParticleSystem('coulomb_crystals/crystal-evolution_' + nCrystal) 
    sol = ODEint(s, np.array([0,0,0]), timeStep, timeStep, ODESystem=ODESystemExact,  Step=StepEulerAdvanced)
    PlotFinalPositions(sol)
    
def SolveParticleSystem():
    initsystem = MakeParticleSystem(0,1)
    #SaveParticleSystem(initsystem, '1')

    #initsystem = LoadParticleSystem('1') 
    #initsystem = LoadParticleSystem('oneElectron')
    
    #trapParams = np.array([0, 0.12, 0.45])
    trapParams = np.array([0, 0.05, 0.6])#test
    
    ODESolution = ODEint(initsystem, trapParams, 1*endTime, timeStep, ODESystemExact, StepEulerAdvanced)
    
    #PlotCoordinates(ODESolution)
    PlotEnergy(ODESolution)
    PlotODESolution(ODESolution)

    
def MakeStabilityDiagram():
    #initsystem = LoadParticleSystem('1_170')
    initsystem = LoadParticleSystem('test')     
    #initsystem = MakeParticleSystem(0,1)
    #SaveParticleSystem(initsystem, 'test')
    
    stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, 0.0, 0.25, 10, 0.0, 1.0, 10, endTime, timeStep) 
    #stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, 0.0, 0.25, 160, 0.0, 1.0, 160, endTime, timeStep) 
        
    PlotStability(stability, plotParams)   

def MakeStabilityDiagramEdge(previousFileName=None, q1Start=0.0, q1Stop=0.04, q1Resol=8, q2Start=0.0, q2Stop=0.5, q2Resol=8):
    initsystem = LoadParticleSystem('test')
    
    start = timer()#to track real time of the computation


    if previousFileName == None:
        stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, endTime, timeStep)
        nParticles = plotParams[-4]
        nIons, nElectrons = nParticles
    else:
        q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nIons, nElectrons, eta = ParseFileName(previousFileName)
        
    for i in range(2):
        
        scaleFactor = 2
        
        previousFileName = str(int(nIons)) + '_ions_' + str(int(nElectrons)) + '_electrons_' + 'q1_' + str(q1Start) + '-' + str(q1Stop) + '_q2_' + str(q2Start) + '-' + str(q2Stop)+ '_' + str(int(q2Resol)) + 'x' + str(int(q1Resol)) + '_' + str(int(f2/f1))

        q1Resol = scaleFactor * q1Resol - (scaleFactor - 1)
        q2Resol = scaleFactor * q2Resol - (scaleFactor - 1)
                
        stability, plotParams = StabilityDiagramEdge(initsystem, ODESystemExact, previousFileName, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, endTime, timeStep)
        
    PlotStability(stability, plotParams)
    
    stop = timer()#to track real time of the computation
    time = stop - start
    
    print('total time: ', time)

    
    
if __name__ == '__main__':
    
    timeStep = 1/(200) #in normal time scale 1/(100 * f2) -> faster period divided into 100 segments 
    endTime = 1*2.5 * f2 / f1 #in normal time scale 5 * (1 / f1) -> five slower periods
    #MakeCoulombCrystal()    
    #TestCoulombCrystal()    
    #PlotCrystalTest() 
    
    #SolveParticleSystem()
    
    #MakeStabilityDiagram()
    MakeStabilityDiagramEdge()    
    
    #params = ClickStabilityRegions('0_ions_1_electrons_q1_0.0-0.04_q2_0.0-0.5_10x10_833')
    #unstableRegion, stableRegion = LoadTriangles(params)

    #SaveTriangles(triangleUnstable, triangleStable, params)
