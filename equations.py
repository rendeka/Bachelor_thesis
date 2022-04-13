import numpy as np
from parameters import * 

"""
stability for light particle is stable for heavy particle if the condition q1 * (f2 / f1)**2 * (m2 / m1) < 0.9
"""


def ODESystemExact(rv, t, aCoulomb, mass, charge, trapParams): #exact equation of motion
    
    #tau = t * f2 / 2
    tau = t
    
    #a = 4 * V0 * charge / (mass * f2**2 * r0**2)
    #q1 = -2 * V1 * charge / (mass * f2**2 * r0**2)
    #q2 = -2 * V2 * charge / (mass * f2**2 * r0**2)
    
    #a=0
    #q1=0
    #q2=0.45
    
    if (mass == electronMass):        
        a, q1, q2 = trapParams
    else:
        a, q1, q2 = trapParams * (electronMass / ionMass)
        
    
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

def ODESystemEffective(rv, t, aCoulomb, mass, charge, trapParams): #effective potential
    
    if (mass == electronMass):        
        a, q1, q2 = trapParams
    else:
        a, q1, q2 = trapParams * (electronMass / ionMass)
    
    r, v = rv
    x,y,z = r
    vx,vy,vz = v

    x1 = vx
    vx1 = aCoulomb[0] - x / 4 * (a + 2 * q1**2 / 2 * (f2 / f1)**2 + q2**2 / 2)
    
    y1 = vy
    vy1 = aCoulomb[1] - y / 4 * (a + 2 * q1**2 / 2 * (f2 / f1)**2 + q2**2 / 2)
    
    z1 = vz
    vz1 = aCoulomb[2] + z / 2 * (a + 2 * q1**2 / 2 * (f2 / f1)**2 + q2**2 / 2)
    #vz1 = aCoulomb[2] - z / 4 * (a + 2 * q1**2 / 2 * (f2 / f1)**2 + q2**2 / 2)#test
    
    r1 = np.array([x1, y1, z1])
    v1 = np.array([vx1, vy1, vz1])

    return np.array([r1, v1])

def ODESystemDampingEff(rv, t, aCoulomb, mass, charge, trapParams): #effective potential
    
    if (mass == electronMass):        
        a, q1, q2 = trapParams
        q1 = 0
    else:
        a, q1, q2 = trapParams * (electronMass / ionMass)
    
    r, v = rv
    x,y,z = r
    vx,vy,vz = v

    x1 = vx
    vx1 = aCoulomb[0] - x / 4 * (a + 2 * q1**2 / 2 * (f2 / f1)**2 + q2**2 / 2) - const2 * vx
    
    y1 = vy
    vy1 = aCoulomb[1] - y / 4 * (a + 2 * q1**2 / 2 * (f2 / f1)**2 + q2**2 / 2) - const2 * vy
    
    z1 = vz
    vz1 = aCoulomb[2] + z / 2 * (a + 2 * q1**2 / 2 * (f2 / f1)**2 + q2**2 / 2) - const2 * vz
    #vz1 = aCoulomb[2] - z / 4 * (a + 2 * q1**2 / 2 * (f2 / f1)**2 + q2**2 / 2) - const2 * vz
    
    r1 = np.array([x1, y1, z1])
    v1 = np.array([vx1, vy1, vz1])

    return np.array([r1, v1])

def ODESystemCrystal(rv, t, aCoulomb, mass, charge, trapParams): 
    
    r, v = rv
    rv1 = [v, aCoulomb - r*const*np.abs(charge)/mass - const2 * v]
    
    return np.array(rv1)

def ODESystemDampingExact(rv, t, aCoulomb, mass, charge, trapParams):

    tau = t
    
    if (mass == electronMass):        
        a, q1, q2 = trapParams
        q1 = 0
    else:
        a, q1, q2 = trapParams * (electronMass / ionMass)
    
    r, v = rv
    x,y,z = r
    vx,vy,vz = v

    x1 = vx
    vx1 = aCoulomb[0] - x * (a - 2 * q1 * np.cos(2 * tau * f1 / f2) - 2 * q2 * np.cos(2 * tau)) - const2 * vx
    
    y1 = vy
    vy1 = aCoulomb[1] - y * (a - 2 * q1 * np.cos(2 * tau * f1 / f2) - 2 * q2 * np.cos(2 * tau)) - const2 * vy
    
    z1 = vz
    vz1 = aCoulomb[2] + 2 * z * (a - 2 * q1 * np.cos(2 * tau * f1 / f2) - 2 * q2 * np.cos(2 * tau)) - const2 * vz
    
    
    r1 = np.array([x1, y1, z1])
    v1 = np.array([vx1, vy1, vz1])

    return np.array([r1, v1])

def GetDt(ODESystem):
    result = 1/200
    return {
        'ODESystemExact': result,
        'ODESystemEffective': result * 10,
        'ODESystemDampingEff': result * 10,
        'ODESystemCrystal': result * 10,
    }[ODESystem.__name__]
