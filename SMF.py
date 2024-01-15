import integration_library as IL
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import classy
from classy import Class
from scipy.integrate import ode, odeint
import scipy.constants as SPC
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
import math
import scipy
from tqdm import tqdm 
import sys
from numpy import diff
import matplotlib.pylab as pl
import matplotlib as mpl
from matplotlib.colors import LogNorm

# Defs

kvec = np.logspace(np.log10(0.00001),np.log10(1000.0),10000)
H0 = 67.66
Omegab0 = 0.02242/(H0/100)**2
Omegam0 = (0.02242/(H0/100)**2+0.11933/(H0/100)**2)
Omegar0 = 8.493e-5
c = 299792.45800000057
Tcmb = 2.72e6
h = 0.6766
c = 3.3
GN = 4.301*10**(-9)
rho = 3*H0**2*Omegam0/(8*np.pi*GN)
rhocr = 2.77536627e11
rhom = rhocr*Omegam0

K0 = [0, 0.5, 1]
beta = np.linspace(0,0.5,10)
H_arr_kmoufl = [[],[],[]]
dH_arr_kmoufl = [[],[],[]]
H_int_kmoufl = [[],[],[]]
dH_int_kmoufl = [[],[],[]]

M_kmoufl = {}


common_settings_kmoufl =  {'n_s':0.9665,
          'A_s':2.101e-9,
          'tau_reio':0.0561,
          'omega_b':0.02242,
          'omega_cdm':0.11933,
          'h':0.6766,
          'YHe':0.2425,
          'T_cmb':2.7255,
          'gauge':'newtonian', #FOR MGCLASS TO WORK, GAUGE NEEDS TO BE NEWTONIAN
          'k_pivot': 0.05,
          'mg_z_init': 111110.000,
          'l_logstep': 1.025,
          'l_linstep':15,
          'P_k_max_1/Mpc':1500.0,
          'l_switch_limber':9,
          'perturb_sampling_stepsize': 0.05,
          'output':'tCl,pCl,lCl,mPk',
          'l_max_scalars': 3000,
          'lensing': 'yes',
          'z_max_pk': 99}

def Pk(a, model, par1, par2):
    common_settings = {'n_s':0.9665,
          'A_s':2.101e-9,
          'tau_reio':0.0561,
          'omega_b':0.02242,
          'omega_cdm':0.11933,
          'h':0.6766,
          'YHe':0.2425,
          'T_cmb':2.7255,
          'gauge':'newtonian', #FOR MGCLASS TO WORK, GAUGE NEEDS TO BE NEWTONIAN
          'k_pivot': 0.05,
          'mg_z_init': 111110.000,
          'l_logstep': 1.025,
          'l_linstep':15,
          'P_k_max_1/Mpc':1500.0,
          'l_switch_limber':9,
          'perturb_sampling_stepsize': 0.05,
          'output':'tCl,pCl,lCl,mPk',
          'l_max_scalars': 3000,
          'lensing': 'yes',
          'z_max_pk': 99}
    common_settings['mg_ansatz'] = model
    if model == 'plk_late':
        common_settings['mg_E11'] = par1
        common_settings['mg_E22'] = par2
    elif model == 'z_flex_late':
        common_settings['mg_muz'] = par1
        common_settings['mg_gamz'] = par2
        common_settings['mg_zzn'] = par2
    elif model == 'z_xpans_late':
        common_settings['mg_T1'] = par1
        common_settings['mg_T2'] = par2
        common_settings['mg_T3'] = par1
        common_settings['mg_T4'] = par2
        common_settings['mg_zzn'] = 1
    elif model == 'GI':
        common_settings['Omega_Lambda'] = 0.0
        common_settings['w0_fld'] = par1
        common_settings['gamGI'] = par2
    elif model == 'nDGP':
        common_settings['rc'] = par1
    else:
        common_settings['beta_kmfl'] = par1
        common_settings['k0_kmfl'] = par2
        
    M = Class()
    M.set(common_settings)
    M.compute()

    Pk = []
    for k in kvec:
        Pk.append(M.pk(k,1/a-1))
    
    return Pk

"""
for i in range(len(K0)):
    for j in tqdm(range(len(beta))):
        common_settings_kmoufl['mg_ansatz'] = 'kmoufl'
        common_settings_kmoufl['beta_kmfl'] = beta[j]
        common_settings_kmoufl['k0_kmfl'] = K0[i]

        M_kmoufl[j] = Class()
        M_kmoufl[j].set(common_settings_kmoufl)
        
        M_kmoufl[j].compute()
        a = np.logspace(-6,0,10000)
        H_arr_kmoufl[i].append([M_kmoufl[j].Hubble(1/ai-1)*c for ai in a])
        dH_arr_kmoufl[i].append(np.gradient(H_arr_kmoufl[i][j])/np.gradient(a))
        H_int_kmoufl[i].append(scipy.interpolate.interp1d(a, H_arr_kmoufl[i][j], fill_value = 'extrapolate'))
        dH_int_kmoufl[i].append(scipy.interpolate.interp1d(a, dH_arr_kmoufl[i][j], fill_value = 'extrapolate'))
"""

def H_f(model_H,a, par1, par2):
    if model_H == 'LCDM':
        return H0*np.sqrt(1-Omegam0-Omegar0+Omegam0*a**(-3)+Omegar0*a**(-4))
    elif model_H == 'wCDM':
        par1 = wL
        return H0*np.sqrt((1-Omegam0-Omegar0)*a**(-3*(1+wL))+Omegam0*a**(-3)+Omegar0*a**(-4))
    elif model_H == 'nDGP':
        rc = par1
        Omegarc = 1/(4*H0**2*rc**2)
        OmegaLambda0 = 1 - Omegam0 - Omegar0 + 2*np.sqrt(Omegarc)
        return H0*np.sqrt(Omegam0*a**(-3)+Omegar0*a**(-4)+Omegarc + OmegaLambda0)-H0*np.sqrt(Omegarc)
    elif model_H == 'kmoufl':
        return H_int_kmoufl[i_kmoufl][j_kmoufl](a)

def dH_f(model_H,a, par1, par2):
    if model_H == 'LCDM':
        return -0.5*(H0*(3*a*Omegam0 + 4*Omegar0))/(a**5*np.sqrt((a*Omegam0 + Omegar0 - a**4*(-1 + Omegam0 + Omegar0))/a**4))
    elif model_H == 'wCDM':
        wL = par1
        return (a**(-5 - 3*wL)*(3*a*H0*(1 + wL)*(-1 + Omegam0 + Omegar0) - a**(3*wL)*H0*(3*a*Omegam0 + 4*Omegar0)))/(2*np.sqrt((Omegar0 + a*(Omegam0 - (-1 + Omegam0 + Omegar0)/a**(3*wL)))/a**4))
    elif model_H == 'nDGP':
        rc = par1
        Omegarc = 1/(4*H0**2*rc**2)
        OmegaLambda0 = 1 - Omegam0 - Omegar0 + 2*np.sqrt(Omegarc)
        return  (H0*((-3*Omegam0)/a**4 - (4*Omegar0)/a**5))/(2.*np.sqrt(OmegaLambda0 + Omegam0/a**3 + Omegar0/a**4 + Omegarc))
    elif model_H == 'kmoufl':
        return dH_int_kmoufl[i_kmoufl][j_kmoufl](a)
        
def mu(model_H, model, par1, par2, a):
    H = H_f(model_H, a, par1, par2)
    dHda = dH_f(model_H, a, par1, par2)
    dHdt = a*H*dHda
    rhom = 3*H0**2*Omegam0*a**(-3)
    rhor = 3*H0**2*Omegar0*a**(-4)
    rhoL = 3*H**2-rhom-rhor
    OmegaL = rhoL/(3*H**2)
    OmegaM = rhom/(3*H**2)
    if model == 'plk_late':
        f1 = par1*OmegaL
        return 1 + f1
    elif model =='z_flex_late':
        return 1+ par1*(1-a)**1 - par1*(1-a)**2
    elif model == 'z_xpans_late':
        return 1 + par1*OmegaL + par2*OmegaL**2
    elif model == 'GI':
        wL = par1
        gamma = par2
        return 2/3*OmegaM**(gamma-1)*(OmegaM**gamma+2-3*gamma+3*(gamma-1/2)*(OmegaM+(1+wL)*OmegaL))
    elif model == 'nDGP':
        rc = par1
        beta = 1 + 2*H/c*par1*(1+dHdt/(3*H**2))
        return 1 + 1/(3*beta)
    elif model == 'kmoufl':
        beta = par1
        K0 = par2
        A_kmfl = 1.0 + par1*a
        X_kmfl = 0.5 * A_kmfl**2*(H*a)**2/((1-Omegam0-Omegar0)*H0**2)
        k_prime_mfl = 1.0 + 2.0*par2*X_kmfl
        epsl1_kmfl = 2.0*par1**2/k_prime_mfl    
        X_kmfl_dot = 0.5 * A_kmfl**2/((1-Omegam0-Omegar0)*H0**2)*2.0*H*a*(H*H*a+dHdt*a)
        return 1.+ epsl1_kmfl
    else:
        sys.exit("Check whether the mu(a) model is correctly chosen")

def delta_nl_ODE(a, y, model_H, model, par1, par2):
    delta,ddeltada = y
    H = H_f(model_H, a, par1, par2)
    dH = dH_f(model_H, a, par1, par2)
    dddeltada          = -(3/a+dH/H)*ddeltada + (3*Omegam0*mu(model_H, model, par1, par2, a))/(2*a**5*H**2/H0**2)*delta*(1+delta) + 4*ddeltada**2/(3*(1+delta))
    return [ddeltada, dddeltada]

def delta_l_ODE(a, y, model_H, model, par1, par2):
    delta,ddeltada = y
    H = H_f(model_H, a, par1, par2)
    dH = dH_f(model_H, a, par1, par2)
    dddeltada          = -(3/a+dH/H)*ddeltada + (3*Omegam0*mu(model_H, model, par1, par2, a))/(2*a**5*H**2/H0**2)*delta
    return [ddeltada, dddeltada]

def collapse(deltai, model_H, model, par1, par2):
    ai = 1e-5
    dt = 0.0001
    ddeltai = deltai/ai
    init = [deltai, ddeltai]
    system = ode(delta_nl_ODE)
    system.set_f_params(*[model_H, model, par1, par2])
    system.set_initial_value(init, ai)
    
    data = [0,0]
    while system.successful() and system.y[0] <= 1e7 and system.t <= 1/(1-0.5):
        data = [system.t + dt, system.integrate(system.t + dt)[0]]   
        #plt.scatter(system.t+dt, system.integrate(system.t + dt)[0], c ='tab:blue')
    data = np.array(data)
    ac = data[0]

    return ac

def interpolate_ac(ac, model_H, model, par1, par2):
    delta_i = np.logspace(-4,-1,1000)
    arr_di_ac = []
    for i in tqdm(range(len(delta_i))):
        arr_di_ac.append([delta_i[i],collapse(delta_i[i], model_H, model,par1,par2)])
    interp_di_ac = scipy.interpolate.interp1d(np.array(arr_di_ac)[:,1],np.array(arr_di_ac)[:,0],fill_value = 'extrapolate')
    return interp_di_ac(ac)

def linear(deltai_collapse, a, model_H, model, par1, par2):
    ai = 1e-5
    dt = 0.0001
    ddeltai = deltai_collapse/ai    
    init = [deltai_collapse, ddeltai]
    system = ode(delta_l_ODE)
    system.set_f_params(*[model_H, model, par1, par2])
    system.set_initial_value(init, ai)
    
    data = []
    while system.successful() and system.t <= a:
        data.append([system.t + dt, system.integrate(system.t + dt)[0]])   
    
    data = np.array(data)
    return data

def delta_c_at_ac(model_H, model, ac, par1, par2):
    return linear(interpolate_ac(ac, model_H, model, par1, par2),ac, model_H, model,par1,par2)[-1,1]


def sigma(k,Pk,R):
    yinit = np.array([0.0], dtype=np.float64)
    eps   = 1e-13  #change this for higher/lower accuracy
    h1    = 1e-12
    hmin  = 0.0
    beta_ST = 4.8
    W   = (1+(k*R)**beta_ST)**(-1)
    Pk1 = Pk*W**2*k**2/(2.0*np.pi**2)
    
    return np.sqrt(IL.odeint(yinit, k[0],k[-1], eps,
                             h1, hmin, np.log10(k), Pk1,
                             'sigma', verbose=False)[0])

def dSdM(k, Pk, rhoM, M):
    c_ST = 3.3
    R1=(3.0*M/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)
    s1=sigma(k,Pk,R1)

    M2=M*1.0001
    R2=(3.0*M2/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)
    s2=sigma(k,Pk,R2)

    return (s2-s1)/(M2-M)


def ST_mass_function(k, Pk, rhoM, Masses, a, model_H, model, par1, par2):
    c_ST = 3.3
    dndM = np.zeros(Masses.shape[0], dtype=np.float64)
    deltac = delta_c_at_ac(model_H, model, a, par1, par2)
    
    for i,M in enumerate(Masses):
        R   = (3.0*M/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)
        nu  = (deltac/sigma(k,Pk,R))**2

        dndM[i]=-(rhoM/M)*dSdM(k,Pk,rhoM,M)/sigma(k,Pk,R)
        dndM[i]*=0.3222*np.sqrt(2*nu/np.pi)*(1+1/(nu**0.3))
        dndM[i]*=np.exp(-0.5*nu)

    return dndM

def func_SFR(x,a):
    z = 1/a-1
    nu = np.exp(-4*a**2)
    alpha_SFR = -1.412+0.731*(a-1)*nu
    Delta_SFR = 3.508 + (2.608*(a-1)-0.043*z)*nu
    gamma_SFR = 0.316+(1.319*(a-1)+0.279*z)*nu
    return -np.log10(10**(alpha_SFR*x)+1)+Delta_SFR*(np.log10(1+np.exp(x)))**gamma_SFR/(1+np.exp(10**(-x)))

def Mstar(Mh, model_SFR, a):
    z = 1/a-1
    if model_SFR == 'phenomenological_regular':
        epstar = 0.32
        Mg = Omegab0/Omegam0*Mh        
        Mstar = epstar*Mg
    elif model_SFR == 'phenomenological_extreme':
        Mg = Omegab0/Omegam0*Mh
        Mstar = Mg
    elif model_SFR == 'Behroozi':
        nu = np.exp(-4*a**2)
        log10M1 = 11.514+(-1.793*(a-1)+(-0.251)*z)*nu
        log10eps = -1.777+(-0.006*(a-1)+(-0.000)*z)*nu-0.119*(a-1)
        log10Mstar = log10eps + log10M1 + func_SFR(np.log10(Mh)-log10M1,a)-func_SFR(0,a)
        Mstar = 10**log10Mstar
    elif model_SFR == 'double_power':
        Mg = Omegab0/Omegam0*Mh
        Mp = 10**12.1
        A = 10**(-1.69)
        gamma = 1.32
        delta = 0.43
        epstar = A*Mh/((Mh/Mp)**(-gamma)+(Mh/Mp)**delta)
        Mstar = epstar*Mh
    else:
        sys.exit("Incorrect SFR model is being used")
    return Mstar

def SMF(Mstar_var,k, Pk, rhoM, a, model_H, model,model_SFR, par1, par2):
    Mh_arr = np.logspace(6.5,16,1000)
    Mstar_arr = Mstar(Mh_arr, model_SFR, a)
    SMF = scipy.interpolate.interp1d(Mstar_arr,ST_mass_function(k,Pk,rhoM,Mh_arr,a,model_H,model,par1,par2)*Mh_arr*np.gradient(np.log10(Mh_arr))/np.gradient(np.log10(Mstar_arr)), fill_value="extrapolate")
    return SMF(Mstar_var)

Masses = np.logspace(6.5,15,50)
model_H = 'LCDM'
model = 'plk_late'
model_SFR = 'Behroozi'
par1 = 0
par2 = 0
a = 1
z = 1/a-1
Masses_star = Mstar(Masses, model_SFR, 1/(1+z)) 
#varepsilon = -1.777+(-0.006*(a-1)+(-0.000)*z)*nu-0.119*(a-1)
plt.plot(Masses_star,h**(-3)*np.gradient(np.log10(Masses_star),np.log10(Masses))*ST_mass_function(kvec/h, np.array(Pk(a, model, par1, par2))*h**3, rhom, Masses, a, model_H, model, par1, par2)*Masses*np.log(10), color = "tab:blue")
#plt.plot(Masses, Masses_star)
x = np.loadtxt('Downloads/z0pt1.dat')[:,0] #https://github.com/bmoster/emerge/blob/master/data/smf.dat
y = np.loadtxt('Downloads/z0pt1.dat')[:,1]
plt.scatter(10**x,10**y)
#y = scipy.interpolate.interp1d(10**x,10**y, fill_value='extrapolate')
#plt.plot(Masses_star,np.gradient(np.log10(Masses_star),np.log10(Masses))*ST_mass_function(kvec/h, np.array(Pk(a, model, par1, par2))*h**3, rhom, Masses, a, model_H, model, par1, par2)*Masses*np.log(10)/y(Masses_star), color = "tab:blue")
#
plt.xscale('log')
plt.yscale('log')
plt.xlim(10**8,10**12)
#plt.ylim(0,6)
plt.savefig('SMF.pdf', bbox_inches = 'tight')
