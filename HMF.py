import integration_library as IL
%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import classy
from classy import Class


common_settings = {'A_s':2.101e-9,
          'n_s':0.9665,
          'tau_reio':0.0561,
          'omega_b':0.02242,
          'omega_cdm':0.11933,
          'h':0.6766,
          'YHe':0.2425,
          'T_cmb':2.7255,
          'gauge':'newtonian', #FOR MGCLASS TO WORK, GAUGE NEEDS TO BE NEWTONIAN
          'k_pivot': 0.05,
          'mg_z_init': 10.000,
          'l_logstep': 1.025,
          'l_linstep':15,
          'P_k_max_1/Mpc':3.0,
          'l_switch_limber':9,
          'perturb_sampling_stepsize': 0.05,
          'output':'tCl,pCl,lCl,mPk',
          'l_max_scalars': 3000,
          'lensing': 'yes',
          'mg_ansatz':'plk_late',
          'mg_E11': 0,
          'mg_E22': 0}

M = Class()
M.set(common_settings)
M.compute()

kvec = np.logspace(np.log10(0.0001),np.log10(1.0),1000)
twopi = 2.*math.pi

Tcmb = 2.72e6
h = 0.6766
c = 3.3
GN = 4.301*10**(-9)
rho = 3*H0**2*Omegam0/(8*np.pi*GN)

Pk = []
for k in kvec:
    Pk.append(M.pk(k,0.))
        
def H_f(a):
    if model_H == 'LCDM':
        return H0*np.sqrt(1-Omegam0-Omegar0+Omegam0*a**(-3)+Omegar0*a**(-4))
    elif model_H == 'wCDM':
        return H0*np.sqrt((1-Omegam0-Omegar0)*a**(-3*(1+wL))+Omegam0*a**(-3)+Omegar0*a**(-4))
    elif model_H == 'nDGP':
        Omegarc = 1/(4*H0**2*rc**2)
        OmegaLambda0 = 1 - Omegam0 - Omegar0 + 2*np.sqrt(Omegarc)
        return H0*np.sqrt(Omegam0*a**(-3)+Omegar0*a**(-4)+Omegarc + OmegaLambda0)-H0*np.sqrt(Omegarc)
    elif model_H == 'kmoufl':
        return H_int_kmoufl[i_kmoufl][j_kmoufl](a)

def dH_f(a):
    if model_H == 'LCDM':
        return -0.5*(H0*(3*a*Omegam0 + 4*Omegar0))/(a**5*np.sqrt((a*Omegam0 + Omegar0 - a**4*(-1 + Omegam0 + Omegar0))/a**4))
    elif model_H == 'wCDM':
        return (a**(-5 - 3*wL)*(3*a*H0*(1 + wL)*(-1 + Omegam0 + Omegar0) - a**(3*wL)*H0*(3*a*Omegam0 + 4*Omegar0)))/(2*np.sqrt((Omegar0 + a*(Omegam0 - (-1 + Omegam0 + Omegar0)/a**(3*wL)))/a**4))
    elif model_H == 'nDGP':
        Omegarc = 1/(4*H0**2*rc**2)
        OmegaLambda0 = 1 - Omegam0 - Omegar0 + 2*np.sqrt(Omegarc)
        return  (H0*((-3*Omegam0)/a**4 - (4*Omegar0)/a**5))/(2.*np.sqrt(OmegaLambda0 + Omegam0/a**3 + Omegar0/a**4 + Omegarc))
    elif model_H == 'kmoufl':
        return dH_int_kmoufl[i_kmoufl][j_kmoufl](a)
        
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
        beta = 1 + 2*H/c*par1*(1+dHdt/(3*H**2))
        return 1 + 1/(3*beta)
    elif model == 'kmoufl':
        A_kmfl = 1.0 + par1*a
        X_kmfl = 0.5 * A_kmfl**2*(H*a)**2/((1-Omegam0-Omegar0)*H0**2)
        k_prime_mfl = 1.0 + 2.0*par2*X_kmfl
        epsl1_kmfl = 2.0*par1**2/k_prime_mfl    
        X_kmfl_dot = 0.5 * A_kmfl**2/((1-Omegam0-Omegar0)*H0**2)*2.0*H*a*(H*H*a+dHdt*a)
        return 1.+ epsl1_kmfl

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


def sigma(k,Pk,R):
    yinit = np.array([0.0], dtype=np.float64)
    eps   = 1e-13  #change this for higher/lower accuracy
    h1    = 1e-12
    hmin  = 0.0
    beta = 4.8
    W   = (1+(k*R)**beta)**(-1)
    Pk1 = Pk*W**2*k**2/(2.0*np.pi**2)
    
    return np.sqrt(IL.odeint(yinit, k[0],k[-1], eps,
                             h1, hmin, np.log10(k), Pk1,
                             'sigma', verbose=False)[0])

def sigma_M(k,Pk,M):
    R =  (M**(1/3)*(3/np.pi)**(1/3))/(2**(2/3)*c*rho**(1/3))
    return np.array(sigma(k, Pk, R))

def dsigma_M(k,Pk,M):
    sigma_arr = np.array([sigma_M(k, Pk, m) for m in M])
    dx = np.gradient(np.log(M))
    dy = np.gradient(sigma_arr)
    return dy/dx

def nu(pars1, pars2, k, Pk, M):
    deltac = linear(find_deltai(1,pars1, pars2),1,pars1, pars2)[-1,1]
    return deltac**2/sigma_M(np.array(k),np.array(Pk),M)**2
    
def fnu(pars1, pars2, k, Pk, M):
    A = 0.3222
    p = 0.3
    nu_v = nu(pars1, pars2, k, Pk, M)
    return np.array(A*np.sqrt(2*nu_v**2/np.pi)*(1+nu_v**(-2*p))*np.exp(-nu_v**2/2))
