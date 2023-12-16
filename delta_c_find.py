from scipy.integrate import ode, odeint
import scipy.constants as SPC
import numpy as np
import matplotlib.pyplot as plt
import math

H0 = 67.66
Omegam0 = (0.02242/(H0/100)**2+0.11933/(H0/100)**2)
Omegar0 = 8.493e-5

def H_f(a):
    if model_H == 'LCDM':
        return H0*np.sqrt(1-Omegam0-Omegar0+Omegam0*a**(-3)+Omegar0*a**(-4))
    elif model_H == 'wCDM':
        return H0*np.sqrt((1-Omegam0-Omegar0)*a**(-3*(1+wL))+Omegam0*a**(-3)+Omegar0*a**(-4))
    elif model_H == 'nDGP':
        Omegarc = 1/(4*H0**2*rc**2)
        OmegaLambda0 = 1 - Omegam0 - Omegar0
        return H0*np.sqrt(Omegam0*a**(-3)+Omegar0*a**(-4)+Omegarc+OmegaLambda0)-H0*np.sqrt(Omegarc)

def dH_f(a):
    if model_H == 'LCDM':
        return -0.5*(H0*(3*a*Omegam0 + 4*Omegar0))/(a**5*np.sqrt((a*Omegam0 + Omegar0 - a**4*(-1 + Omegam0 + Omegar0))/a**4))
    elif model_H == 'wCDM':
        return (a**(-5 - 3*wL)*(3*a*H0*(1 + wL)*(-1 + Omegam0 + Omegar0) - a**(3*wL)*H0*(3*a*Omegam0 + 4*Omegar0)))/(2*np.sqrt((Omegar0 + a*(Omegam0 - (-1 + Omegam0 + Omegar0)/a**(3*wL)))/a**4))
    elif model_H == 'nDGP':
        Omegarc = 1/(4*H0**2*rc**2)
        OmegaLambda0 = 1 - Omegam0 - Omegar0
        return -0.5*(H0*(3*a*Omegam0 + 4*Omegar0))/(a**5*np.sqrt((a*Omegam0 + Omegar0)/a**4 + Omegarc + OmegaLambda0))
        
def mu(par1, par2, a):
    H = H_f(a)
    dHda = dH_f(a)
    dHdt = a*H*dHda
    rhom = 3*H0**2*Omegam0*a**(-3)
    rhor = 3*H0**2*Omegar0*a**(-4)
    rhoL = 3*H**2-rhom-rhor
    OmegaL = rhoL/(3*H**2)
    OmegaM = rhom/(3*H**2)
    if model == 'E11':
        f1 = par1*OmegaL
        return 1 + f1
    elif model =='gmu':
        return 1+ par1*(1-a)**1 - par1*(1-a)**2
    elif model == 'DES':
        return 1 + par1*OmegaL + par2*OmegaL**2
    elif model == 'wCDM':
        return 2/3*OmegaM**(gamma-1)*(OmegaM**gamma+2-3*gamma+3*(gamma-1/2)*(OmegaM+(1+wL)*OmegaL))
    elif model == 'nDGP':
        OmegaLambda0 = 1 - Omegam0 - Omegar0
        Omegarc = 1/(4*H0**2*par1**2)
        beta = 1 + (Omegam0*a**(-3)+Omegar0*a**(-4)+2*OmegaLambda0)/(2*np.sqrt(Omegarc*(Omegam0*a**(-3)+Omegar0*a**(-4)+OmegaLambda0)))
        return 1 + 1/(3*beta)
        
def delta_nl_ODE(a, y, par1, par2):
    delta,ddeltada = y
    H = H_f(a)
    dH = dH_f(a)
    dddeltada          = -(3/a+dH/H)*ddeltada + (3*Omegam0*mu(par1, par2, a))/(2*a**5*H**2/H0**2)*delta*(1+delta) + 4*ddeltada**2/(3*(1+delta))
    return [ddeltada, dddeltada]

def delta_l_ODE(a, y, par1, par2):
    delta,ddeltada = y
    H = H_f(a)
    dH = dH_f(a)
    dddeltada          = -(3/a+dH/H)*ddeltada + (3*Omegam0*mu(par1, par2, a))/(2*a**5*H**2/H0**2)*delta
    return [ddeltada, dddeltada]

def collapse(deltai, par1, par2):
    ai = 1e-5
    dt = 0.01
    ddeltai = deltai/ai
    init = [deltai, ddeltai]
    system = ode(delta_nl_ODE)
    system.set_f_params(*[par1, par2])
    system.set_initial_value(init, ai)
    
    data = [0,0]
    while system.successful() and system.y[0] <= 1e7 and system.t <= 1:
        data = [system.t + dt, system.integrate(system.t + dt)[0]]   
        #plt.scatter(system.t+dt, system.integrate(system.t + dt)[0], c ='tab:blue')
    data = np.array(data)
    ac = data[0]

    return ac

def find_deltai(ac_input, par1, par2):
    deltarange = np.logspace(-5,-2,num=1000)
    deltai_collapse = deltarange[-1]
    for i in range(len(deltarange)):
        a_v = collapse(deltarange[i], par1, par2)
        b_v = ac_input
        if abs(a_v - b_v)/((a_v+b_v)/2) <= 0.01:
            deltai_collapse = deltarange[i]
    return deltai_collapse

def linear(deltai_collapse, a, par1, par2):
    ai = 1e-5
    dt = 0.01
    ddeltai = deltai_collapse/ai    
    init = [deltai_collapse, ddeltai]
    system = ode(delta_l_ODE)
    system.set_f_params(*[par1, par2])
    system.set_initial_value(init, ai)
    
    data = []
    while system.successful() and system.t <= a:
        data.append([system.t + dt, system.integrate(system.t + dt)[0]])   
        #plt.scatter(system.t+dt, system.integrate(system.t + dt)[0], c ='tab:orange')
    
    data = np.array(data)
    return data
