import numpy as np

def StepEuler(ODESystem, rv, t, dt, aCoulomb, mass, charge, trapParams):
    
    rv1 = rv + ODESystem(rv, t, aCoulomb, mass, charge, trapParams) * dt
    t1 = t + dt
    
    return rv1, t1

def StepEuler2(ODESystem, rv, t, dt, aCoulomb, mass, charge, trapParams):
    
    k1 = ODESystem(rv, t, aCoulomb, mass, charge, trapParams)
    k2 = ODESystem(rv + k1 * dt, t + dt, aCoulomb, mass, charge, trapParams)
    
    rv1 = rv + 0.5 * (k1 + k2) * dt
    t1 = t + dt
    
    return rv1, t1

def StepEulerAdvanced(ODESystem, rv, t, dt, aCoulomb, mass, charge, trapParams):
    
    k1 = ODESystem(rv, t, aCoulomb, mass, charge, trapParams)
    k2 = ODESystem(rv + k1 * 0.5 * dt, t + 0.5 * dt, aCoulomb, mass, charge, trapParams)
    
    rv1 = rv + k2 * dt
    t1 = t + dt
    
    return rv1, t1

def StepVerlet(ODESystem, rv, t, dt, aCoulomb, mass, charge, trapParams):

    r, v = rv
    v, a = ODESystem(rv, t, aCoulomb, mass, charge, trapParams)
    
    r1 = r + v * dt + 0.5 * a * dt**2
    
    a1 = ODESystem(np.array([r1, v]), t, aCoulomb, mass, charge, trapParams)[1]    
    
    v1 = v + 0.5 * (a + a1) * dt
    t1 = t + dt
    
    rv1 = np.array([r1, v1])
    
    return rv1, t1   

def StepRungaKutta(ODESystem, rv, t, dt, aCoulomb, mass, charge, trapParams):
    
    k1 = ODESystem(rv, t, aCoulomb, mass, charge, trapParams)
    k2 = ODESystem(rv + 0.5 * k1 * dt, t + 0.5 * dt, aCoulomb, mass, charge, trapParams)
    k3 = ODESystem(rv + 0.5 * k2 * dt, t + 0.5 * dt, aCoulomb, mass, charge, trapParams)
    k4 = ODESystem(rv + k3 * dt, t + dt, aCoulomb, mass, charge, trapParams)
    
    rv1 = rv + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4) * dt
    t1 = t + dt 
    
    return rv1, t1