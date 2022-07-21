"""This file is concerned with ploting and saving pictures"""
import numpy as np
import matplotlib.pyplot as plt


from matplotlib.cm import ScalarMappable
from matplotlib import cm # Colour maps for the contour graph
from matplotlib import rc # Latex typesetting

from scipy import interpolate

from fileHandling import *

from mpl_toolkits import mplot3d
from parameters import *

#rc('font', **{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True) 


def PlotODESolution(ODESolution):
    """plots trajectories of all particles in the system"""    
    rs, _, methodName, exeTime, _, system, _ = ODESolution
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    n = len(system)
    
    for i in range(n):
            
        x = rs[i,:,0] * 1000 # 1000 to convert into milimeters
        y = rs[i,:,1] * 1000
        z = rs[i,:,2] * 1000
        
        charge = system[i][3]
        if charge > 0:
            color = "red"
        else:
            color = "blue"
        
        #plt.axis('off')
        #plt.grid(b=None)
        ax.plot3D(x, y, z, color)
    
    ax.view_init(15, 240)
    ax.set_title(str(n) + " charged particles \n" + "method: " + methodName + "\n" + "execution time: " + str(exeTime) + " s")
    
    ax.set_xlabel('x[mm]')
    ax.set_ylabel('y[mm]')
    ax.set_zlabel('z[mm]')
    
    ax.xaxis.set_tick_params(labelsize=sizeTick)
    ax.yaxis.set_tick_params(labelsize=sizeTick)
    ax.zaxis.set_tick_params(labelsize=sizeTick)
    
    framesize = r0 * 500
    
    ax.auto_scale_xyz([-framesize, framesize], [-framesize, framesize], [-framesize, framesize])


    plt.grid(b=None)
    #plt.legend()
    
    fileName = 'energies'
    extensions = ['eps', 'png']
    
    for extension in extensions: #saving pictures 
        plt.savefig('pics/crystal/' + fileName + str(methodName) + '.' + extension, dpi=500, format=extension)

    plt.show()   
    
def PlotODESolution2D(ODESolution):
    """plots trajectories of all particles in the system"""    
    rs, _, methodName, exeTime, _, system, _ = ODESolution
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    n = len(system)
    
    for i in range(n):
            
        x = rs[i,:,0] * 1000 # 1000 to convert into milimeters
        y = rs[i,:,1] * 1000
        
        charge = system[i][3]
        if charge > 0:
            color = "red"
        else:
            color = "blue"
        
        #plt.axis('off')
        #plt.grid(b=None)
        ax.plot(x, y, color)
    
    #ax.view_init(15, 240)
    ax.set_title(str(n) + " charged particles \n" + "method: " + methodName + "\n" + "execution time: " + str(exeTime) + " s")
    
    ax.set_xlabel('x[mm]')
    ax.set_ylabel('y[mm]')
    
    ax.xaxis.set_tick_params(labelsize=sizeTick)
    ax.yaxis.set_tick_params(labelsize=sizeTick)
    
    #plt.gca().set_aspect('equal') #setting aspect ration to 1

    plt.grid(b=None)
    #plt.legend()
    plt.savefig(fname='pics/' + str(methodName), dpi=500)
    plt.show()   

def PlotEnergy(ODESolution):
    """plots energy balance of all particles throughout the computation"""
    
    _, _, methodName, _, energies, _, _ = ODESolution
    
    totalEnergy = energies[0]
    
    fig = plt.figure()

    #plt.title("Total energy of the system in time")
    plt.xlabel('t')
    plt.ylabel('Energy')    
        
    plt.plot(np.arange(len(totalEnergy)),energies[0] ,label='Total energy')
    plt.plot(np.arange(len(totalEnergy)),energies[1] ,label='Kinetic energy')
    plt.plot(np.arange(len(totalEnergy)),energies[2] ,label='Potential energy')
    
    plt.legend()
    
    fileName = 'energies'
    extensions = ['eps', 'png']
    
    for extension in extensions: #saving pictures 
        plt.savefig('pics/crystal/' + fileName + '.' + extension, format=extension)

    plt.show()
    
def PlotFinalPositions(ODESolution):
    """scatter plot of all particle positions -> using this for coulomb crystals """
    
    rs, vs, methodName, exeTime, energy, system, _ = ODESolution
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    n = len(system)
    
    for i in range(n):
            
        """
        x = rs[i,-1,0] * 1000 # 1000 to convert into milimeters
        y = rs[i,-1,1] * 1000 
        z = rs[i,-1,2] * 1000 
        """
        
        x = rs[i,0] * 1000 # 1000 to convert into milimeters
        y = rs[i,1] * 1000 
        z = rs[i,2] * 1000 
        
        charge = system[-1][3]
        if charge > 0:
            col = "red"
        else:
            col = "blue"

        ax.scatter(x, y, z, color=col)
    
    ax.view_init(15, 240)
    
    ax.set_xlabel('x[mm]', fontsize=sizeLabelSmall)
    ax.set_ylabel('y[mm]', fontsize=sizeLabelSmall)
    ax.set_zlabel('z[mm]', fontsize=sizeLabelSmall)
    
    ax.xaxis.set_tick_params(labelsize=sizeTickSmall)
    ax.yaxis.set_tick_params(labelsize=sizeTickSmall)
    ax.zaxis.set_tick_params(labelsize=sizeTickSmall)
    
    framesize = 0.5
    ax.auto_scale_xyz([-framesize, framesize], [-framesize, framesize], [-framesize, framesize])
    
    
    plt.grid(b=None)
    #plt.legend()
    #plt.savefig('pics/' + "final-positions", dpi=500)
    fileName = 'final-positions'
    extensions = ['eps', 'png']
    
    for extension in extensions: #saving pictures 
        plt.savefig('pics/crystal/' + fileName + '.' + extension, dpi=500, format=extension)
    plt.show()
    
def PlotStability(data=np.zeros((2,2)), params=np.zeros(10), index=None, fileName=None, velocityDiagram=False):
    
    if fileName is None:
        q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, time, f1, f2 = params
        fileName = MakeFileName(params)

    else:
        data, params = LoadStabilityDiagram(fileName=fileName, velocityDiagram=velocityDiagram)
        q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, eta, f1, f2 = params
                
    x_vals = np.linspace(q2Start, q2Stop, int(q2Resol))
    y_vals = np.linspace(q1Start, q1Stop, int(q1Resol))
    X, Y = np.meshgrid(x_vals, y_vals) 
            
    if index != None:
        fileName = fileName + '_' + str(index)
            
    path = 'pics/stability_diagrams/'
    if velocityDiagram:
        path = path + 'velocity/'
        fig, ax = plt.subplots()
        
        vmin, vmax = 0, 50
        levels = 400 
        level_boundaries = np.linspace(vmin, vmax, levels + 1)
        
        quadcontourset = ax.contourf(
            X, Y, data,
            levels,
            vmin=vmin, vmax=vmax
        )
        
        
        cbar = fig.colorbar(
            ScalarMappable(norm=quadcontourset.norm, cmap=quadcontourset.cmap),
            ticks=range(vmin, vmax+5, 15),
            boundaries=level_boundaries,
            values=(level_boundaries[:-1] + level_boundaries[1:]) / 2,
            extend='max',
        )
        
        cbar.ax.tick_params(labelsize=sizeTick)
        cbar.set_label(r'$\dfrac{\bar{\mathcal{v}}}{\mathcal{v}_0}$', fontsize=sizeLabel, rotation=0)
        
        
    else:
        fig = plt.figure()        
        plt.contourf(X, Y, data)

    
    plt.xlabel('$q_{2}$')
    plt.ylabel('$q_{1}$')
    
    if velocityDiagram:
        extensions = ['png', 'pdf']
    else:    
        extensions = ['eps', 'png', 'pdf']    
    for extension in extensions: #saving pictures 
        plt.savefig(path + fileName + '.' + extension, format=extension)
    
    plt.show()
    
def PlotStabilityRescaled(data=np.zeros((2,2)), params=[], index=None, fileName=None, velocityDiagram=False):
        
    if fileName is None:
        q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, f1, f2 = params
        fileName = MakeFileName(params)

    else:
        data, params = LoadStabilityDiagram(fileName, velocityDiagram=velocityDiagram)
        #data = RepairData(data)
        if fileName[0] == 'd':            
            q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, eta, f1, f2 = params
        else:
            q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, eta, f1, f2 = params            
                
    x_vals = np.linspace(q2Start, q2Stop, int(q2Resol)) * 10 
    y_vals = np.linspace(q1Start, q1Stop, int(q1Resol)) * 100
    X, Y = np.meshgrid(x_vals, y_vals) 
            
    if index != None:
        fileName = fileName + '_' + str(index)
            
    path = 'pics/stability_diagrams/'
    if velocityDiagram:
        #path = path + 'velocity/'
        fig, ax = plt.subplots()
        
        vmin, vmax = 0, 21
        levels = 1000    
        level_boundaries = np.linspace(vmin, vmax, levels + 1)
        
        quadcontourset = ax.contourf(
            X, Y, data,
            levels,
            vmin=vmin, vmax=vmax
        )
        
        
        cbar = fig.colorbar(
            ScalarMappable(norm=quadcontourset.norm, cmap=quadcontourset.cmap),
            ticks=range(vmin, vmax+5, 5),
            boundaries=level_boundaries,
            values=(level_boundaries[:-1] + level_boundaries[1:]) / 2,
            extend='max',
        )
        
        cbar.ax.tick_params(labelsize=sizeTickSmall)
        cbar.set_label(r'$\dfrac{\langle v \rangle}{v_0}$', fontsize=sizeLabelSmall, rotation=0, labelpad=15)
        
        
    else:
        fig = plt.figure()        
        plt.contourf(X, Y, data)
        
        
    plt.xticks(fontsize=sizeTickSmall)    
    plt.yticks(fontsize=sizeTickSmall)    
    plt.xlabel('$q_{2} \ [10^{-1}]$', fontsize=sizeLabelSmall)
    plt.ylabel('$q_{1} \ [10^{-2}]$', fontsize=sizeLabelSmall)
    plt.tight_layout()
    
    if velocityDiagram:
        extensions = ['png', 'pdf']
    else:    
        extensions = ['eps', 'png', 'pdf']
    
    for extension in extensions: #saving pictures 
        plt.savefig(path + 'rescaled/' + fileName + '.' + extension, format=extension)
    
    plt.show()

"""    
def PlotStabilityVelocity(data=np.zeros((2,2)), params=np.zeros(10), index=None, fileName=None):
    
    if fileName is None:
        q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, time, f1, f2 = params
    else:
        data, params = LoadStabilityDiagram(fileName)
        q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, eta, f1, f2 = params
        
    
    x_vals = np.linspace(q2Start, q2Stop, int(q2Resol))
    y_vals = np.linspace(q1Start, q1Stop, int(q1Resol))
    X, Y = np.meshgrid(x_vals, y_vals)
    
    fig, ax = plt.subplots()
    
    #vmin, vmax = 0, 40
    vmin, vmax = 0, 100
    levels = 1000
    level_boundaries = np.linspace(vmin, vmax, levels + 1)
    
    quadcontourset = ax.contourf(
        X, Y, data,
        levels,
        vmin=vmin, vmax=vmax
    )
    
    
    fig.colorbar(
        ScalarMappable(norm=quadcontourset.norm, cmap=quadcontourset.cmap),
        ticks=range(vmin, vmax+5, 5),
        boundaries=level_boundaries,
        values=(level_boundaries[:-1] + level_boundaries[1:]) / 2,
        extend='max',
    )
      
    plt.xlabel('$q_{2}$')
    plt.ylabel('$q_{1}$')
    
"""
"""
    fileName = MakeFileName(params)
    
    if index != None:
        fileName = fileName + '_' + str(index)
        
    extensions = ['eps', 'png']
    
    for extension in extensions: #saving pictures 
        plt.savefig('pics/stability_diagrams/' + fileName + '.' + extension, format=extension)
    
    plt.show()    
"""

def PlotCrystalTest(nCrystal='20'):
    s = LoadParticleSystem('coulomb_crystals/crystal-evolution_' + nCrystal) 
    sol = ODEint(s, np.array([0,0,0]), 1e-10, 1e-10, ODESystem=ODESystemExact,  Step=StepEulerAdvanced, freezeIons=True)
    PlotFinalPositions(sol)
    
def PlotCrystal(nCrystal='20'):
    s = LoadParticleSystem('coulomb_crystals/' + nCrystal) 
    #sol = ODEint(s, np.array([0,0,0]), 1e-10, 1e-10, ODESystem=ODESystemExact,  Step=StepEulerAdvanced, freezeIons=True)
    rs=[]
    for p in s:
        rs.append(p[0])
    rs = np.array(rs)
    sol = rs, None, None, None, None, s, None
    PlotFinalPositions(sol)   
    
    
def PlotStabilityDet(data, params):
    
    q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, eta = params
                
    x_vals = np.linspace(q2Start, q2Stop, int(q2Resol))
    y_vals = np.linspace(q1Start, q1Stop, int(q1Resol))
    X, Y = np.meshgrid(x_vals, y_vals)
        
    fileName = 'q1_' + str(q1Start) + '-' + str(q1Stop) + '_q2_' + str(q2Start) + '-' + str(q2Stop) + '_' + str(int(q2Resol)) + 'x' + str(int(q1Resol)) + '_' + str(eta)
 
    path = 'pics/stability_diagrams/determinant/'
    
    fig = plt.figure()        
    plt.contourf(X, Y, data)
    
    plt.xlabel('$q_{2}$')
    plt.ylabel('$q_{1}$')
    
    extensions = ['eps', 'png', 'pdf']
    
    for extension in extensions: #saving pictures 
        plt.savefig(path + fileName + '.' + extension, format=extension)
    
    plt.show()
    
def RepairData(data, defect=1.0, fix=velocityChop):
    for i, row in enumerate(data):
        for j, element in enumerate(row):
            #if (i != 0)or(j != 0):
            if element == defect:
                data[i,j] = fix
    return data