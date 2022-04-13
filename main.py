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
    s=LoadParticleSystem('coulomb_crystals/crystal-evolution_' + nCrystal) 
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
    initsystem = LoadParticleSystem('1')     
    #initsystem = MakeParticleSystem(0,1)
    #SaveParticleSystem(initsystem, '1_170')
    stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, 0, 0.14, 20, 0, 0.6, 20, endTime, timeStep) 

    PlotStability(stability, plotParams)    

if __name__ == '__main__':
    
    timeStep = 1/(200) #in normal time scale 1/(100 * f2) -> faster period divided into 100 segments 
    endTime = 2*2.5 * f2 / f1 #in normal time scale 5 * (1 / f1) -> five slower periods

    #MakeCoulombCrystal()    
    #TestCoulombCrystal()    
    #PlotCrystalTest() 
    
    #SolveParticleSystem()
    
    #MakeStabilityDiagram()    
    
    #"""
    triangleStable = []
    triangleUnstable = []

    data, params = LoadStabilityDiagram('0_ions_1_electrons_q1_0-0.14_q2_0-0.6_20x20_13')
    
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
    #"""

