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
    
def SolveParticleSystem(nCrystal='20'):
    """doesn't work if number of electrons=0 and freezedIons=True"""    
    
    #initsystem = MakeParticleSystem(50,1)
    #initsystem = MakeParticleSystem(0,1)
    
    #initsystem = np.array([], dtype=object)    
    #initsystem = LoadParticleSystem('test')
    
    #initsystem = LoadParticleSystem('coulomb_crystals/crystal-evolution_' + nCrystal)
    #initsystem = AddElectron(initsystem)     

    #initsystem = LoadParticleSystem('1') 
    #initsystem = LoadParticleSystem('oneElectron')
    
    
    #"""
    crystal = LoadParticleSystem('coulomb_crystals/' + nCrystal) 
    electron = LoadParticleSystem('oneElectron')[0]
    #electron = MakeParticle(electronMass, electronCharge, 100)
    initsystem = AddElectron(crystal, electron)  
    initsystem = np.array([electron])
    #SaveParticleSystem(initsystem, 'oneElectron')
    #"""
    
    #trapParams = np.array([0, 0.24, 0.37])
    trapParams = np.array([0, 0.02, 0.35])
    
    ODESolution = ODEint(initsystem, trapParams, 1*endTime, 1*timeStep, ODESystemExact, StepVerlet, freezeIons=True)
    
    PlotEnergy(ODESolution)
    PlotODESolution2D(ODESolution)
    PlotODESolution(ODESolution)

    
def MakeStabilityDiagram(nCrystal='20', q1Start=0.0, q1Stop=0.1, q1Resol=9, q2Start=0.0, q2Stop=0.55, q2Resol=8, freezeIons=True):
    
    #initsystem = LoadParticleSystem('1_170')
    
    initsystem = LoadParticleSystem('wtf') 
    
    #initsystem = MakeParticleSystem(0,1)
    #SaveParticleSystem(initsystem, 'wtf')
    
    #electron = LoadParticleSystem('oneElectron')[0]
    
    """
    initsystem = np.array([electron], dtype=object)
    stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, endTime, timeStep, freezeIons)         
    PlotStability(stability, plotParams) 
    """
    
    """
    #crystal = LoadParticleSystem('coulomb_crystals/crystal' + '50') 
    crystal = LoadParticleSystem('coulomb_crystals/50') 
    electron = LoadParticleSystem('oneElectron')[0]
    initsystem = AddElectron(crystal, electron)
    SaveParticleSystem(initsystem, 'test')
    """

    
    """
    initsystem = MakeCoulombCrystal(nCrystal='5', trapParams=np.array([0, 0.02, 0.35]))
    initsystem = AddElectron(initsystem, electron)
    SaveParticleSystem(initsystem, 'test')
    """
    
    freezeIons = True
    stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, endTime, timeStep, freezeIons)       
    PlotStabilityVelocity(stability, plotParams)
    #PlotStability(stability, plotParams)   
    
def MakeStabilityDiagramList(nCrystal='20', q1Start=0.0, q1Stop=0.05, q1Resol=12, q2Start=0.0, q2Stop=0.5, q2Resol=12, freezeIons=True):
    
    systems = CreateInitialSystems()
    
    freezeIons = True

    for idx, initsystem in enumerate(systems):    
        stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, endTime, timeStep, freezeIons)       
        PlotStability(stability, plotParams, idx)   
    
    
def MakeStabilityDiagramEdge(nCrystal='20', previousFileName=None, q1Start=0.0, q1Stop=0.1, q1Resol=12, q2Start=0.0, q2Stop=0.55, q2Resol=12, freezeIons=True):
    
    initsystem = LoadParticleSystem('wtf')
    #initsystem = MakeParticleSystem(0,1)
    #SaveParticleSystem(initsystem, 'wtf')

    #initsystem = LoadParticleSystem('coulomb_crystals/crystal-evolution_' + nCrystal)
    #initsystem = AddElectron(initsystem)
    
    #electron = LoadParticleSystem('oneElectron')[0]
    #initsystem = np.array([electron], dtype=object)
    
    #crystal = LoadParticleSystem('coulomb_crystals/' + nCrystal) 
    #electron = LoadParticleSystem('oneElectron')[0]
    #initsystem = np.array([electron])

    #initsystem = AddElectron(crystal, electron)
    #SaveParticleSystem(initsystem, '50stability')
    
    start = timer()#to track real time of the computation
    freezeIons = True

    if previousFileName is None:
        stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, endTime, timeStep, freezeIons)
        nParticles = plotParams[-4]
        nIons, nElectrons = nParticles
    else:
        q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nIons, nElectrons, eta = ParseFileName(previousFileName)
        
    for _ in range(6):
        print('resolution: ' +  str(int(q1Resol)) + 'x' + str(int(q2Resol)))
        for i in range(2):#after this loop resolution will be doubled in both coordinates q1 adn q2
            subStart = timer()
        
            params = q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, (nIons, nElectrons), int(f2/f1), f1, f2            
            fileName = MakeFileName(params)
            
            if i == 0: scaleFactor = 3/2 #we want to avoid scaling resolution by factor of 2, because in that case new points
            if i == 1: scaleFactor = 4/3 #in the grid would be exactly between the old ones making unnecessarily big rounding errors
            
            q1Resol = int(scaleFactor * q1Resol)# we must cast into integer in order to allow any racional number as scaleFactor
            q2Resol = int(scaleFactor * q2Resol)# but it's your responsibility to make sure that resolution has no decimal places for every interation
                    
            stability, plotParams = StabilityDiagramEdge(initsystem, ODESystemExact, fileName, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, endTime, timeStep, freezeIons)            
            PlotStability(stability, plotParams)
            
            
            print('time: ', timer() - subStart)
    
    stop = timer()#to track real time of the computation
    time = stop - start
    
    print('total time: ', time)
    
def StabilityDiagramForCrystals(q1Start=0.0, q1Stop=0.1, q1Resol=32, q2Start=0.0, q2Stop=1.0, q2Resol=32):
    
    crystal = LoadParticleSystem('coulomb_crystals/crystal-evolution_' + '20') 
    electron = MakeParticle(electronMass, electronCharge)
    initsystem = np.array([], dtype=object)
    initsystem = AddElectron(initsystem, electron)
    
    for idx in range(1):
                
        freezeIons = True
        stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, endTime, timeStep, freezeIons)           
        PlotStability(stability, plotParams)   
        
        initsystem = AddElectron(crystal, electron)        
        freezeIons = True
        stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, endTime, timeStep, freezeIons)           
        PlotStability(stability, plotParams) 
        #print('testing result ', idx)    
        #TestCoulombCrystal(nCrystal=str(idx), trapParams=tp)

def PlotStabilityEdge(fileName, plotParams=None):
    fMatrix, _, _ = PrepForStabilityEdge(fileName)
    params = ParseFileName(fileName)
    q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nIons, nElectrons, eta = params
    
    for iX, fX in enumerate(fMatrix):
        for iY, fY in enumerate(fX):
            if fY == -1:
                fMatrix[iX,iY] = 1
            else:
                fMatrix[iX,iY] = -1
        
    x_vals = np.linspace(q2Start, q2Stop, int(q2Resol))
    y_vals = np.linspace(q1Start, q1Stop, int(q1Resol))
    X, Y = np.meshgrid(x_vals, y_vals)
        
    vmin = 0
    vmax = 1000
    
    levels = np.linspace(vmin, vmax, 2+1)
    fig,ax = plt.subplots()
    contourf_ = ax.contourf(X,Y,fMatrix, levels=levels, vmin=vmin, vmax=vmax)
    #cbar = fig.colorbar(contourf_)
    
    #plt.show()
    
def PlotSVelocityEdge(fileName1, fileName2):
    fMatrix, _, _ = PrepForStabilityEdge(fileName2)
    data, params = LoadStabilityDiagram(fileName1)
    q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, eta, f1, f2 = params
    
    for iX, fX in enumerate(fMatrix):
        for iY, fY in enumerate(fX):
            if fY == -1:
                fMatrix[iX,iY] = 1
            else:
                fMatrix[iX,iY] = -1
        
    x_vals = np.linspace(q2Start, q2Stop, int(q2Resol))
    y_vals = np.linspace(q1Start, q1Stop, int(q1Resol))
    X, Y = np.meshgrid(x_vals, y_vals)
    
    fig, (ax1, ax2) = plt.subplots(nrows=2)
    
    vmin1, vmax1 = 0, 100
    levels1= 1000
    level_boundaries = np.linspace(vmin1, vmax1, levels1 + 1)
    
    quadcontourset = ax1.contourf(
        X, Y, data,
        levels1,
        vmin=vmin1, vmax=vmax1
    )
    
    
    fig.colorbar(
        ScalarMappable(norm=quadcontourset.norm, cmap=quadcontourset.cmap),
        ticks=range(vmin1, vmax1+5, 5),
        boundaries=level_boundaries,
        values=(level_boundaries[:-1] + level_boundaries[1:]) / 2,
    )
      
    vmin2 = 0
    vmax2 = 1000
    
    levels2 = np.linspace(vmin2, vmax2, 2+1)
    contourf_ = ax2.contourf(X,Y,fMatrix, levels=levels2, vmin=vmin2, vmax=vmax2)
    #cbar = fig.colorbar(contourf_)
    
    plt.xlabel('$q_{2}$')
    plt.ylabel('$q_{1}$')
    
    plt.show()
    
def ThePlot(fileName1, fileName2):
    PlotStabilityVelocity(fileName=fileName1)
    PlotStabilityEdge(fileName=fileName2)
    plt.show()
        
 
if __name__ == '__main__':
    
    timeStep = 1/(200) #in normal time scale 1/(100 * f2) -> faster period divided into 100 segments 
    endTime = 2*2.5 * f2 / f1 #in normal time scale 5 * (1 / f1) -> five slower periods
    
    tp = np.array([0, 0.02, 0.35])
    #MakeCoulombCrystal(nCrystal='20', trapParams=tp)
    #print('testing result')    
    #TestCoulombCrystal(nCrystal='20', trapParams=tp)   
    
    #MakeCoulombCrystalFromGrid(nCrystal='19', trapParams=tp)
    
    #PlotCrystalTest(nCrystal='20') 
    #PlotCrystal(nCrystal='5') 
    
    #SolveParticleSystem(nCrystal='20')
    
    #MakeStabilityDiagram()  
    
    #PlotStabilityEdge(fileName='0_ions_1_electrons_q1_0.0-0.1_q2_0.0-0.55_768x897_3')
    #PlotSVelocityEdge(fileName='0_ions_1_electrons_q1_0.0-0.05_q2_0.0-0.5_512x512_170')
    #PlotStabilityVelocity(fileName='0_ions_1_electrons_q1_0.0-0.1_q2_0.0-0.55_8x9_3')
    
    PlotSVelocityEdge(fileName1='0_ions_1_electrons_q1_0.0-0.1_q2_0.0-0.55_192x192_3', fileName2='0_ions_1_electrons_q1_0.0-0.1_q2_0.0-0.55_768x897_3')

    #MakeStabilityDiagramList()
        
    #MakeStabilityDiagramEdge(nCrystal='50', previousFileName='50_ions_1_electrons_q1_0.0-0.05_q2_0.0-0.5_128x128_13')   
    #MakeStabilityDiagramEdge()   
    
    #StabilityDiagramForCrystals()
    
    #params = ClickStabilityRegions('0_ions_1_electrons_q1_0.0-0.04_q2_0.0-0.46_12x12_170')
    #unstableRegion, stableRegion = LoadTriangles(params)

    #SaveTriangles(triangleUnstable, triangleStable, params)
