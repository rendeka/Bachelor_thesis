import numpy as np
from parameters import * 

"""
This file contains systems of ODEs that we want to integrate
"""


def ODESystemExact(rv, tau, aCoulomb, mass, charge, trapParams): #exact equation of motion
    """Exact equations of motion for ideal quadrupole trap""" 
    
    a, q1, q2 = trapParams * (electronMass / mass)#trap parameters depend on charge to mass ration  
        
    r, v = rv
    x,y,z = r
    vx,vy,vz = v
            
    x1 = vx
    vx1 = aCoulomb[0] - x * (a - 2 * q1 * np.cos(2 * tau * f1 / f2) - 2 * q2 * np.cos(2 * tau))
    
    y1 = vy
    vy1 = aCoulomb[1] - y * (a - 2 * q1 * np.cos(2 * tau * f1 / f2) - 2 * q2 * np.cos(2 * tau))
    
    z1 = vz
    vz1 = aCoulomb[2] + 2 * z * (a - 2 * q1 * np.cos(2 * tau * f1 / f2) - 2 * q2 * np.cos(2 * tau))
    #vz1 = aCoulomb[2] - z * (a - 2 * q1 * np.cos(2 * tau * f1 / f2) - 2 * q2 * np.cos(2 * tau))#test
    
    r1 = np.array([x1, y1, z1])
    v1 = np.array([vx1, vy1, vz1])

    return np.array([r1, v1])

def ODESystemExactSymmetric(rv, tau, aCoulomb, mass, charge, trapParams):
    """Symmetrized equations of motion for ideal quadrupole trap""" 
    
    a, q1, q2 = trapParams * (electronMass / mass)#trap parameters depend on charge to mass ration  

    r, v = rv
            
    r1 = v
    v1 = aCoulomb - r * (a - 2 * q1 * np.cos(2 * tau * f1 / f2) - 2 * q2 * np.cos(2 * tau))

    return np.array([r1, v1])

def ODESystemEffective(rv, t, aCoulomb, mass, charge, trapParams):
    """Equations of motion devived from effective potential for ideal quadrupole trap""" 
 
    if (mass == electronMass):        
        a, q1, q2 = trapParams
        q1 = 0
    else:
        a, q1, q2 = trapParams * (electronMass / mass)    
    r, v = rv
    x,y,z = r
    vx,vy,vz = v

    x1 = vx
    vx1 = aCoulomb[0] - x / 4 * (a + q1**2 / 2 * (f2 / f1)**2 + q2**2 / 2)
    
    y1 = vy
    vy1 = aCoulomb[1] - y / 4 * (a + q1**2 / 2 * (f2 / f1)**2 + q2**2 / 2)
    
    z1 = vz
    vz1 = aCoulomb[2] - z * (a + q1**2 / 2 * (f2 / f1)**2 + q2**2 / 2)
    
    r1 = np.array([x1, y1, z1])
    v1 = np.array([vx1, vy1, vz1])

    return np.array([r1, v1])

def ODESystemEffectiveSymmetric(rv, tau, aCoulomb, mass, charge, trapParams): 
    """Symmetrized equations of motion devived from effective potential for ideal quadrupole trap""" 
    
    if (mass == electronMass):        
        a, q1, q2 = trapParams
        q1 = 0
    else:
        a, q1, q2 = trapParams * (electronMass / mass)
        
    r, v = rv

    r1 = v
    v1 = aCoulomb - r / 4 * (a + 2 * q1**2 / 2 * (f2 / f1)**2 + q2**2 / 2)

    return np.array([r1, v1])

def ODESystemEffectiveDamping(rv, tau, aCoulomb, mass, charge, trapParams): #effective potential
    """Equations of motion devived from effective potential for ideal quadrupole trap with damping """ 

    if (mass == electronMass):        
        a, q1, q2 = trapParams
        q1 = 0
    else:
        a, q1, q2 = trapParams * (electronMass / mass)
    
    r, v = rv
    x,y,z = r
    vx,vy,vz = v

    x1 = vx
    vx1 = aCoulomb[0] - x / 4 * (a + q1**2 / 2 * (f2 / f1)**2 + q2**2 / 2) - 2 * beta * vx
    
    y1 = vy
    vy1 = aCoulomb[1] - y / 4 * (a + q1**2 / 2 * (f2 / f1)**2 + q2**2 / 2) - 2 * beta * vy
    
    z1 = vz
    vz1 = aCoulomb[2] - z * (a + q1**2 / 2 * (f2 / f1)**2 + q2**2 / 2) - 2 * beta * vz
    
    r1 = np.array([x1, y1, z1])
    v1 = np.array([vx1, vy1, vz1])

    return np.array([r1, v1])

def ODESystemDampingExact(rv, tau, aCoulomb, mass, charge, trapParams):
    
    a, q1, q2 = trapParams * (electronMass / mass)#trap parameters depend on charge to mass ration  
    
    r, v = rv
    x,y,z = r
    vx,vy,vz = v

    x1 = vx
    vx1 = aCoulomb[0] - x * (a - 2 * q1 * np.cos(2 * tau * f1 / f2) - 2 * q2 * np.cos(2 * tau)) - 2 * beta * vx
    
    y1 = vy
    vy1 = aCoulomb[1] - y * (a - 2 * q1 * np.cos(2 * tau * f1 / f2) - 2 * q2 * np.cos(2 * tau)) - 2 * beta * vy
    
    z1 = vz
    vz1 = aCoulomb[2] + 2 * z * (a - 2 * q1 * np.cos(2 * tau * f1 / f2) - 2 * q2 * np.cos(2 * tau)) - 2 * beta * vz
    
    
    r1 = np.array([x1, y1, z1])
    v1 = np.array([vx1, vy1, vz1])

    return np.array([r1, v1])

def ODESystemCrystal(rv, t, aCoulomb, mass, charge, trapParams): 
    
    r, v = rv
    rv1 = [v, aCoulomb - r*const*np.abs(charge)/mass]
    
    return np.array(rv1)

def GetDt(ODESystem):
    result = 1/200
    return {
        'ODESystemExact': result,
        'ODESystemExactSymmetric': result,
        'ODESystemDampingExact': result,
        'ODESystemEffective': result * 10,
        'ODESystemEffectiveDamping': result * 10,
        'ODESystemCrystal': result * 10,
    }[ODESystem.__name__]
