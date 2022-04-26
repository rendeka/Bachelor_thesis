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
from fileHandling import *
from stability import *

def MakeCoulombCrystal(nCrystal=20, trapParams=np.array([0, 0.4, 0])):
    """
    CLUSTER
    Makes coulomb crystal using molecular dynamics. Its single argument final number of ions in crystal
    """
    system, solution = IonCrystal(nIons=1, trapParams=trapParams)
    
    while len(system) < nCrystal:
        
        #input("Press Enter to continue...")
        print('Number of ions: ', len(system))
        system = AddIon(system)
        system, solution = IonCrystal(system=system, trapParams=trapParams)
        
    system, solution = IonCrystal(system=system, trapParams=trapParams)
    
    SaveParticleSystem(system, 'coulomb_crystals/' + str(int(nCrystal)))
        
    #PlotODESolution(solution)
    #PlotFinalPositions(solution)
    
    
def TestCoulombCrystal(nCrystal='20', trapParams=np.array([0, 0.4, 0])):
    system = LoadParticleSystem('coulomb_crystals/' + nCrystal)
    trapParams = np.array([0, 0.4, 0])
    
    for particle in system:
        particle[1] = np.zeros(3)
        
    solution = ODEint(system, trapParams, 1*endTime, timeStep, ODESystem=ODESystemExact,  Step=StepEulerAdvanced)      

    print('total velocity of the system', TotalVelocity(system))
    
    SaveParticleSystem(system, 'coulomb_crystals/crystal-evolution_' + nCrystal)
    
def PlotCrystalTest(nCrystal='20'):
    s = LoadParticleSystem('coulomb_crystals/crystal-evolution_' + nCrystal) 
    sol = ODEint(s, np.array([0,0,0]), timeStep, timeStep, ODESystem=ODESystemExact,  Step=StepEulerAdvanced)
    PlotFinalPositions(sol)
    
def SolveParticleSystem():
    #initsystem = MakeParticleSystem(0,1)
    #SaveParticleSystem(initsystem, '1')
    #initsystem = LoadParticleSystem('test')

    nCrystal = '20'    
    initsystem = LoadParticleSystem('coulomb_crystals/crystal-evolution_' + nCrystal)
    initsystem = AddElectron(initsystem)     

    #initsystem = LoadParticleSystem('1') 
    #initsystem = LoadParticleSystem('oneElectron')
    
    #trapParams = np.array([0, 0.12, 0.45])
    trapParams = np.array([0, 0.06, 0.4])
    
    ODESolution = ODEint(initsystem, trapParams, 3*endTime, timeStep, ODESystemExact, StepEulerAdvanced)
    
    #PlotCoordinates(ODESolution)
    #PlotEnergy(ODESolution)
    PlotODESolution(ODESolution)

    
def MakeStabilityDiagram():
    
    #initsystem = LoadParticleSystem('1_170')
    
    #initsystem = LoadParticleSystem('test') 
    
    #initsystem = MakeParticleSystem(0,1)
    #SaveParticleSystem(initsystem, 'wtf')
    
    nCrystal = '20'    
    initsystem = LoadParticleSystem('coulomb_crystals/crystal-evolution_' + nCrystal)
    initsystem = AddElectron(initsystem)
    
    stability, plotParams = StabilityDiagram(initsystem, ODESystemExactSymmetric, 0.0, 0.14, 20, 0.0, 0.9, 20, endTime, timeStep) 
        
    PlotStability(stability, plotParams)   

def MakeStabilityDiagramEdge(previousFileName=None, q1Start=0.0, q1Stop=0.14, q1Resol=128, q2Start=0.0, q2Stop=0.9, q2Resol=128):
    #initsystem = LoadParticleSystem('test')
    #initsystem = MakeParticleSystem(0,1)
    #SaveParticleSystem(initsystem, 'wtf')

    nCrystal = '20'    
    initsystem = LoadParticleSystem('coulomb_crystals/crystal-evolution_' + nCrystal)
    initsystem = AddElectron(initsystem)
    
    start = timer()#to track real time of the computation

    if previousFileName == None:
        stability, plotParams = StabilityDiagram(initsystem, ODESystemExactSymmetric, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, endTime, timeStep)
        nParticles = plotParams[-4]
        nIons, nElectrons = nParticles
    else:
        q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nIons, nElectrons, eta = ParseFileName(previousFileName)
        
    for _ in range(3):
        for i in range(2):
            params = q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, (nIons, nElectrons), int(f2/f1), f1, f2            
            #previousFileName = str(int(nIons)) + '_ions_' + str(int(nElectrons)) + '_electrons_' + 'q1_' + str(q1Start) + '-' + str(q1Stop) + '_q2_' + str(q2Start) + '-' + str(q2Stop)+ '_' + str(int(q2Resol)) + 'x' + str(int(q1Resol)) + '_' + str(int(f2/f1))
            previousFileName = MakeFileName(params)
            
            if i == 0: scaleFactor = 3/2
            if i == 1: scaleFactor = 4/3
            
            q1Resol = int(scaleFactor * q1Resol)# we must cast into integer in order to allow any racional number as scaleFactor
            q2Resol = int(scaleFactor * q2Resol)# but it's your responsibility to make sure that resolution has no decimal places for every interation
                    
            stability, plotParams = StabilityDiagramEdge(initsystem, ODESystemExactSymmetric, previousFileName, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, endTime, timeStep)
            
            PlotStability(stability, plotParams)
    
    stop = timer()#to track real time of the computation
    time = stop - start
    
    print('total time: ', time)
   
if __name__ == '__main__':
    
    timeStep = 1/(200) #in normal time scale 1/(100 * f2) -> faster period divided into 100 segments 
    endTime = 2*2.5 * f2 / f1 #in normal time scale 5 * (1 / f1) -> five slower periods
    
    #tp = np.array([0, 0.02, 0.38])
    #MakeCoulombCrystal(trapParams=tp)
    #print('testing result')    
    #TestCoulombCrystal(trapParams=tp)    
    
    #PlotCrystalTest() 
    
    SolveParticleSystem()
    
    #MakeStabilityDiagram()
    #MakeStabilityDiagramEdge()   
    
    #params = ClickStabilityRegions('0_ions_1_electrons_q1_0.0-0.04_q2_0.0-0.5_768x768_833')
    #unstableRegion, stableRegion = LoadTriangles(params)

    #SaveTriangles(triangleUnstable, triangleStable, params)
