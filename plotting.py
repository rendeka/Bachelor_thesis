import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm           # Colour maps for the contour graph
from scipy import interpolate

from mpl_toolkits import mplot3d
from parameters import * 


def PlotODESolution(ODESolution):
    
    rs, _, methodName, exeTime, _, system, _ = ODESolution
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    n = len(system)
    
    for i in range(n):
            
        x = rs[i,:,0] * 1000#* 2 / f2 * 1000
        y = rs[i,:,1] * 1000 #* 2 / f2 * 1000
        z = rs[i,:,2] * 1000 #* 2 / f2 * 1000
        
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
    
    ax.xaxis.set_tick_params(labelsize=10)
    ax.yaxis.set_tick_params(labelsize=10)
    ax.zaxis.set_tick_params(labelsize=10)

    plt.grid(b=None)
    #plt.legend()
    plt.savefig(fname='pics/' + str(methodName), dpi=500)
    plt.show()    

def PlotEnergy(ODESolution):
    
    _, _, methodName, _, energies, _, _ = ODESolution
    
    totalEnergy = energies[0]

    #plt.title("Total energy of the system in time")
    plt.xlabel('t')
    plt.ylabel('Energy')    
        
    plt.plot(np.arange(len(totalEnergy)),energies[0] ,label='Total energy')
    plt.plot(np.arange(len(totalEnergy)),energies[1] ,label='Kinetic energy')
    plt.plot(np.arange(len(totalEnergy)),energies[2] ,label='Potential energy')
    
    plt.legend()

    plt.show()
    
def PlotCoordinates(ODESolution):
    
    rs, vs, methodName, exeTime, energy, _, _ = ODESolution
    coordinate = ["x", "y", "z"]
    
    for i in range(2):  
        r = rs[0,:,i] * 2 / f2
        
        plt.title("jednotlive suradnice")
        plt.xlabel = 't'
        plt.ylabel = coordinate[i]            
        plt.plot(np.arange(len(energy)),r ,label=methodName)
        plt.show()
    
def PlotFinalPositions(ODESolution):
    
    rs, vs, methodName, exeTime, energy, system, _ = ODESolution
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    n = len(system)
    
    for i in range(n):
            
        x = rs[i,-1,0] * 1000 #* 2 / f2 * 1000
        y = rs[i,-1,1] * 1000 #* 2 / f2 * 1000
        z = rs[i,-1,2] * 1000 #* 2 / f2 * 1000
        
        charge = system[-1][3]
        if charge > 0:
            color = "red"
        else:
            color = "blue"

        ax.scatter(x, y, z, color)
    
    ax.view_init(15, 240)
    
    ax.set_xlabel('x[mm]')
    ax.set_ylabel('y[mm]')
    ax.set_zlabel('z[mm]')
    
    ax.xaxis.set_tick_params(labelsize=10)
    ax.yaxis.set_tick_params(labelsize=10)
    ax.zaxis.set_tick_params(labelsize=10)

    plt.grid(b=None)
    #plt.legend()
    #plt.savefig('pics/' + "final-positions", dpi=500)
    plt.show()
    
def PlotStability(data, params):
    
    q1Start, q1Stop, q1Resol, q2Start, q2Stop, q2Resol, nParticles, time, f1, f2 = params
    
    plt.figure()
        
    x_vals = np.linspace(q2Start, q2Stop, int(q2Resol))
    y_vals = np.linspace(q1Start, q1Stop, int(q1Resol))
    X, Y = np.meshgrid(x_vals, y_vals) 
    
    
    #cp = plt.contourf(X, Y, data)#, cmap=cm.viridis)
    #plt.colorbar(cp)
    plt.xlabel('$q_{2}$')
    plt.ylabel('$q_{1}$')
    #plt.ylabel("Calorie Burnage")
    #plt.ylabel = '$q_{1}$'
    #plt.clabel('$q_{2}$')
    plt.contourf(X, Y, data)
    
    nIons, nElectrons = nParticles
    fileName = str(int(nIons)) + '_ions_' + str(int(nElectrons)) + '_electrons_' + 'q1_' + str(q1Start) + '-' + str(q1Stop) + '_q2_' + str(q2Start) + '-' + str(q2Stop)+ '_' + str(int(q2Resol)) + 'x' + str(int(q1Resol)) + '_' + str(int(f2/f1))
  
    extensions = ['eps', 'png']
    
    #plt.style.use('seaborn-colorblind')

    
    for extension in extensions:      
        plt.savefig('pics/stability_diagrams/' + fileName + '.' + extension, format=extension)
    
    plt.show()