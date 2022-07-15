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
from time import localtime, strftime

import warnings

    
def SolveParticleSystem(nCrystal='30'):
    """doesn't work if number of electrons=0 and freezedIons=True"""    
    
    #initsystem = MakeParticleSystem(50,1)
    n = int(nCrystal)
    #initsystem = MakeParticleSystem(n,0)
    
    initsystem = LoadParticleSystem('this_3')
    
    #initsystem = LoadParticleSystem('coulomb_crystals/crystal-evolution_' + nCrystal)
    #initsystem = AddElectron(initsystem)     

    #initsystem = LoadParticleSystem('1') 
    #initsystem = LoadParticleSystem('oneElectron')
    
    
    """
    crystal = LoadParticleSystem('coulomb_crystals/' + nCrystal) 
    electron = LoadParticleSystem('oneElectron')[0]
    #electron = MakeParticle(electronMass, electronCharge, 100)
    initsystem = AddElectron(crystal, electron)  
    initsystem = np.array([electron])
    #SaveParticleSystem(initsystem, 'oneElectron')
    """
    
    trapParams = np.array([0, 0.01, 0.25])
    #trapParams = np.array([0, 0.0245, 0.447])
    
    ODESolution = ODEint(initsystem, trapParams, endTime, timeStep, ODESystemExact, StepEulerAdvanced, freezeIons=False)
    #ODESolution = ODEint(initsystem, trapParams, 50*endTime, 5000*timeStep, ODESystemEffectiveDamping, StepEulerAdvanced, freezeIons=False)
    newsystem = ODESolution[-2]
    
    #SaveParticleSystem(newsystem, 'coulomb_crystals/' + str(len(newsystem)))
    PlotODESolution(ODESolution)
    PlotFinalPositions(ODESolution)
    PlotEnergy(ODESolution)

    
def MakeStabilityDiagram(nCrystal='20', q1Start=0.0, q1Stop=0.1, q1Resol=300, q2Start=0.0, q2Stop=0.5, q2Resol=300,
                         freezeIons=True, velocityDiagram=False):
    
    #initsystem = LoadParticleSystem('1_170')
    
    initsystem = LoadParticleSystem('this_3') 
    
    #initsystem = MakeParticleSystem(0,1)
    
    #electron = LoadParticleSystem('oneElectron')[0]
    
    """
    initsystem = np.array([electron], dtype=object)
    stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, q1Start, q1Stop, q1Resol, 
                                             q2Start, q2Stop, q2Resol, endTime, timeStep, freezeIons)         
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
    
    stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, endTime, timeStep, freezeIons, velocityDiagram)       
    PlotStability(data=stability, params=plotParams, velocityDiagram=velocityDiagram)
    
def MakeStabilityDiagramList(nCrystal='20', q1Start=0.0, q1Stop=0.05, q1Resol=40, 
                             q2Start=0.0, q2Stop=0.5, q2Resol=40, freezeIons=True):
    
    systems = CreateInitialSystems()
    
    freezeIons = True

    for idx, initsystem in enumerate(systems):    
        stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, endTime, timeStep, freezeIons)       
        PlotStability(stability, plotParams, idx)   
    
    
def MakeStabilityDiagramEdge(nCrystal='20', previousFileName=None, q1Start=0.0, q1Stop=0.1, q1Resol=80, 
                             q2Start=0.0, q2Stop=0.5, q2Resol=80, freezeIons=True):
    
    #initsystem = LoadParticleSystem('slowElectron')
    #initsystem = LoadParticleSystem('thisOne')
    
    #initsystem = MakeParticleSystem(0,1)
    #SaveParticleSystem(initsystem, 'newSystem')
    
    initsystem = LoadParticleSystem('this_3')
    #initsystem = LoadParticleSystem('large_init_pos_2')
    #initsystem = LoadParticleSystem('large_init_vel')


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

    if previousFileName is None:
        stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, q1Start, q1Stop, q1Resol, 
                                                 q2Start, q2Stop, q2Resol, endTime, timeStep, freezeIons)
        nParticles = plotParams[-4]
        nIons, nElectrons = nParticles
    else:
        q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nIons, nElectrons, eta = ParseFileName(previousFileName)
        
    for _ in range(3):#number of times we will scale the resolution
        print('resolution: ' +  str(int(q1Resol)) + 'x' + str(int(q2Resol)))
        for i in range(1):#after this loop resolution will be doubled in both coordinates q1 adn q2
            subStart = timer()
        
            params = q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, (nIons, nElectrons), int(f2/f1), f1, f2            
            previousFileName = MakeFileName(params)
            
            scaleFactor = 2
            
            q1Resol = int(scaleFactor * q1Resol)# we want to cast into integer in order to allow any racional number as scaleFactor
            q2Resol = int(scaleFactor * q2Resol)# but it's your responsibility to make sure that resolution has no decimal places for every interation
                    
            stability, plotParams = StabilityDiagramEdge(initsystem, ODESystemExact, previousFileName, q1Start, q1Stop, q1Resol, 
                                                         q2Start, q2Stop, q2Resol, endTime, timeStep, freezeIons)            
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
        stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, q1Start, q1Stop, q1Resol,
                                                 q2Start, q2Stop, q2Resol, endTime, timeStep, freezeIons)           
        PlotStability(stability, plotParams)   
        
        initsystem = AddElectron(crystal, electron)        
        freezeIons = True
        stability, plotParams = StabilityDiagram(initsystem, ODESystemExact, q1Start, q1Stop, q1Resol, 
                                                 q2Start, q2Stop, q2Resol, endTime, timeStep, freezeIons)           
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
    
    plt.show()
    
def PlotVelocityEdge(fileNameVel, fileNameStability):
    
    fMatrix, _, _ = PrepForStabilityEdge(fileNameStability)
    q1ResolF, q2ResolF = np.shape(fMatrix)
    
    data, params = LoadStabilityDiagram(fileNameVel, velocityDiagram=True)
    q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, eta, f1, f2 = params
    
    reshapeParams = [q1Start, q1Stop, q1Resol, q1ResolF, q2Start, q2Stop, q2Resol, q2ResolF]
    reshapedData = ReshapeData(data, reshapeParams)
    
    for iX, fX in enumerate(fMatrix):
        for iY, fY in enumerate(fX):
            if fY == -1:
                fMatrix[iX,iY] = 1
            else:
                fMatrix[iX,iY] = -1
        
    x_vals = np.linspace(q2Start, q2Stop, int(q2ResolF)) * 10
    y_vals = np.linspace(q1Start, q1Stop, int(q1ResolF)) * 100
    X, Y = np.meshgrid(x_vals, y_vals)
    
    fig, ax = plt.subplots()
    
    vmin1, vmax1 = 0, 21
    levels1= 1000
    level_boundaries = np.linspace(vmin1, vmax1, levels1 + 1)
    
    quadcontourset = ax.contourf(
        X, Y, reshapedData,
        levels1,
        vmin=vmin1, vmax=vmax1
    )
    
    
    cbar = fig.colorbar(
        ScalarMappable(norm=quadcontourset.norm, cmap=quadcontourset.cmap),
        ticks=range(vmin1, vmax1+5, 5),
        boundaries=level_boundaries,
        values=(level_boundaries[:-1] + level_boundaries[1:]) / 2,
        extend='max',
    )
    
    cbar.ax.tick_params(labelsize=sizeTickSmall)
    cbar.set_label(r'$\dfrac{\langle v \rangle}{v_0}$', fontsize=sizeLabelSmall, rotation=0, labelpad=15)
      
    vmin2 = 0
    vmax2 = 10000
    
    levels2 = np.linspace(vmin2, vmax2, 2+1)
    ax.contourf(X,Y,fMatrix, levels=levels2, vmin=vmin2, vmax=vmax2, cmap=plt.get_cmap('spring'))
    #contourf_ = ax2.contourf(X,Y,fMatrix, levels=levels2, vmin=vmin2, vmax=vmax2)
    #cbar = fig.colorbar(contourf_)
    
    plt.xticks(fontsize=sizeTickSmall)    
    plt.yticks(fontsize=sizeTickSmall)    
    plt.xlabel('$q_{2} \ [10^{-1}]$', fontsize=sizeLabelSmall)
    plt.ylabel('$q_{1} \ [10^{-2}]$', fontsize=sizeLabelSmall)
    plt.tight_layout()
    
    
    fileName = fileNameStability + '_' + str(levels1)
    path = 'pics/stability_diagrams/velocity_with_edge/'
    extensions = ['png', 'pdf']    
    for extension in extensions: #saving pictures 
        plt.savefig(path + fileName + '.' + extension, format=extension)
    
    plt.show()        

if __name__ == '__main__':
    
    warnings.filterwarnings('ignore') # might be a good idea to remove this line when debbuging
    
    startMain = timer()
    
    print(strftime("%b-%d-%Y %H:%M:%S", localtime()))
    
    timeStep = 1/(20)#in normal time scale 1/(10 * f2) -> faster period divided into five segments 
    endTime = 5 * f2 / f1#in normal time scale 50 * (1 / f1) -> fifty slower periods
    
    
    
    #tp = np.array([0, 0.02, 0.35])
    #MakeCoulombCrystal(nCrystal='1', trapParams=tp)
    #print('testing result')    
    #TestCoulombCrystal(nCrystal='20', trapParams=tp)   
    
    #MakeCoulombCrystalFromGrid(nCrystal='30', trapParams=tp, tmax=25*f2/f1*5, dt=1/(20)*5000)
    
    #PlotCrystalTest(nCrystal='20') 
    #PlotCrystal(nCrystal='5') 
    
    #SolveParticleSystem(nCrystal='50')
    
    #PlotStability(fileName='0_ions_1_electrons_q1_0.0-0.05_q2_0.0-0.5_40x41_3',velocityDiagram=True)
    
    #PlotStabilityEdge(fileName='0_ions_1_electrons_q1_0.0-0.1_q2_0.0-0.55_768x897_3')
    #PlotStabilityVelocity(fileName='0_ions_1_electrons_q1_0.0-0.1_q2_0.0-0.5_200x200_3')
    
    #PlotVelocityEdge(fileNameVel='0_ions_1_electrons_q1_0.0-0.05_q2_0.0-0.5_340x340_13', 
    #                    fileNameStability='0_ions_1_electrons_q1_0.0-0.05_q2_0.0-0.5_960x960_13')
    
    #PlotStabilityRescaled(fileName='0_ions_1_electrons_q1_0.01396-0.01626_q2_0.24294-0.24798_240x240_833', velocityDiagram=False)

    #MakeStabilityDiagramList()
        
    #MakeStabilityDiagramEdge(nCrystal='50', previousFileName='50_ions_1_electrons_q1_0.0-0.05_q2_0.0-0.5_128x128_13')   
    MakeStabilityDiagramEdge(previousFileName=None, q1Start=0.0147, q1Stop=0.0156, q1Resol=120, q2Start=0.243, q2Stop=0.249, q2Resol=120)   
    #MakeStabilityDiagram(q1Start=0.0, q1Stop=0.1, q1Resol=360, q2Start=0.0, q2Stop=0.5, q2Resol=360, velocityDiagram=True)  

    #StabilityDiagramForCrystals()
    
    #params = ClickStabilityRegions(fileName='0_ions_1_electrons_q1_0.0-0.05_q2_0.0-0.5_960x960_833')
    """left-click unstable, right-click stable. After you choose triangles you must save them with SaveTriangles"""
    #SaveTriangles(triangleUnstable, triangleStable, params)
    
    #StabilityDet(q1Start=0.0, q1Stop=0.07, q1Resol=990, q2Start=0.0, q2Stop=0.5, q2Resol=990, f1=f1, f2=f2)
    stopMain = timer()
    print(round(stopMain - startMain,2),'seconds')

