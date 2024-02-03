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
from multiprocessing import Pool
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
import scipy.integrate as intg
from scipy.interpolate import InterpolatedUnivariateSpline as _spline

# Defs

kvec = np.logspace(np.log10(0.00001),np.log10(1000.0),10000)
H0 = 67.66
Omegam0 = (0.02242/(H0/100)**2+0.11933/(H0/100)**2)
Omegab0 = 0.02242/(H0/100)**2
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
          'mg_z_init': 10.000,
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
          'mg_z_init': 10.000,
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
def sigma_M(k,Pk,rhoM,M):
    c_ST = 3.3
    R=(3.0*M/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)
    return sigma(k,Pk,R)

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
    deltac = delta_c_at_ac(model_H, model, a, par1, par2)

    if hasattr(Masses, '__len__') and (not isinstance(Masses, str)):
        dndM = np.zeros(Masses.shape[0], dtype=np.float64)
        for i,M in enumerate(Masses):
            R   = (3.0*M/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)
            nu  = (deltac/sigma(k,Pk,R))**2
    
            dndM[i]=-(rhoM/M)*dSdM(k,Pk,rhoM,M)/sigma(k,Pk,R)
            dndM[i]*=0.3222*np.sqrt(2*nu/np.pi)*(1+1/(nu**0.3))
            dndM[i]*=np.exp(-0.5*nu)
    else:
        R   = (3.0*Masses/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)
        nu  = (deltac/sigma(k,Pk,R))**2

        dndM=-(rhoM/Masses)*dSdM(k,Pk,rhoM,Masses)/sigma(k,Pk,R)
        dndM*=0.3222*np.sqrt(2*nu/np.pi)*(1+1/(nu**0.3))
        dndM*=np.exp(-0.5*nu)


    return dndM
    
    
def func_SFR(x,a):
    z = 1/a-1
    nu = np.exp(-4*a**2)
    alpha_SFR = -1.412+0.731*(a-1)*nu
    Delta_SFR = 3.508 + (2.608*(a-1)-0.043*z)*nu
    gamma_SFR = 0.316+(1.319*(a-1)+0.279*z)*nu
    return -np.log10(10**(alpha_SFR*x)+1)+Delta_SFR*(np.log10(1+np.exp(x)))**gamma_SFR/(1+np.exp(10**(-x)))

def epsilon(Mh, model_SFR, a, f0):
    z = 1/a-1
    if model_SFR == 'phenomenological_regular':
        if z<10:
            epstar = 0.15 - 0.03*(z-6)
        else:
            epstar = 0.03
    elif model_SFR == 'phenomenological_extreme':
        epstar = 1
    elif model_SFR == 'Behroozi':
        nu = np.exp(-4*a**2)
        log10M1 = 11.514+(-1.793*(a-1)+(-0.251)*z)*nu
        log10eps = -1.777+(-0.006*(a-1)+(-0.000)*z)*nu-0.119*(a-1)
        log10Mstar = log10eps + log10M1 + func_SFR(np.log10(Mh)-log10M1,a)-func_SFR(0,a)
        Mstar = 10**log10Mstar
        epstar = (Mstar/Mh)/(Omegab0/Omegam0)
    elif model_SFR == 'double_power':
        Mp = 2.8*10**11
        alo = 0.49
        ahi = -0.61
        epstar = 2*f0/((Mh/Mp)**alo + (Mh/Mp)**ahi)
    else:
        sys.exit("Incorrect SFR model is being used")
    return epstar

def varepsilon(Mh, model_SFR, a):
    z = 1/a-1
    if model_SFR == 'phenomenological_extreme' or model_SFR == 'phenomenological_regular':
        varepsilon = 1
    elif model_SFR == 'Behroozi':
        nu = np.exp(-4*a**2)
        log10M1 = 11.514+(-1.793*(a-1)+(-0.251)*z)*nu
        log10eps = -1.777+(-0.006*(a-1)+(-0.000)*z)*nu-0.119*(a-1)
        log10Mstar = log10eps + log10M1 + func_SFR(np.log10(Mh)-log10M1,a)-func_SFR(0,a)
        alpha_SFR = -1.412+0.731*(a-1)*nu
        Delta_SFR = 3.508 + (2.608*(a-1)-0.043*z)*nu
        gamma_SFR = 0.316+(1.319*(a-1)+0.279*z)*nu
        varepsilon = -((Mh**alpha_SFR*alpha_SFR)/(10**(log10M1*alpha_SFR) + Mh**alpha_SFR))+ (np.log(1 + Mh**(1/np.log(10))/np.e**log10M1)**(-1 + gamma_SFR)*Delta_SFR*(10**log10M1*np.exp(10**log10M1/Mh)*(np.e**log10M1 + Mh**(1/np.log(10)))*np.log(10)*np.log(1 + Mh**(1/np.log(10))/np.e**log10M1) + (1 + np.exp(10**log10M1/Mh))*Mh**(1 + 1/np.log(10))*gamma_SFR))/((1 + np.exp(10**log10M1/Mh))**2*Mh*(np.e**log10M1 + Mh**(1/np.log(10)))*np.log(10)**gamma_SFR) 
    elif model_SFR == 'double_power':
        z_interpolate = np.linspace(0,10,11)
        A_int = np.array([-1.69,-1.72,-1.72,-1.27,-1.61,-1.11,-0.8,-1.2,-0.92,-1.34,-1.63])
        delta = 0.43
        Mc = 10**12.3
        gamma_int = np.array([1.32,1.45,1.4,0.61,0.89,1.16,0.88,1.25,1.12,1.52,1.02])
        A = scipy.interpolate.interp1d(z_interpolate,A_int,fill_value='extrapolate')
        gamma = scipy.interpolate.interp1d(z_interpolate,gamma_int,fill_value='extrapolate')
        A = 10**A(z)
        gamma = gamma(z)
        varepsilon = 1- ((-gamma*(Mh/Mc)**(-gamma)+delta*(Mh/Mc)**(delta))/((Mh/Mc)**(-gamma)+(Mh/Mc)**(delta)))
    return varepsilon

def sigma_P(z):
    sigma0 = 0.07
    sigmaz = 0.05
    return sigma0 + sigmaz*z

def f_passive_obs(Masses_star, a):
    z = 1/a-14.86881834
    return ((Masses_star/(10**(10.2+0.5*z)))**(-1.3)+1)**(-1)
    
def SMF_obs(Masses, rhoM, a, model_H, model,model_SFR, par1, par2):
    #Mh_arr = np.logspace(6.5,18,1000)
    HMF = np.log(10)*Masses/h*ST_mass_function(kvec/h, np.array(Pk(a, model, par1, par2))*h**3, rhom, Masses/h, a, model_H, model, par1, par2)/h**3
    Masses_star = epsilon(Masses, model_SFR, a, f0)*Omegab0/Omegab0*Masses
    SMF = HMF*np.gradient(np.log10(Masses))/np.gradient(np.log10(Masses_star)) #varepsilon(Mh_arr,model_SFR,a)

    """
    mu_SMF = -0.020+0.081*(a-1)
    kappa_SMF = 0.045 + (-0.155)*(a-1)
    ci = 0.273*(1+np.exp(1.077-z))**(-1)
    ci1 = 0.273*(1+np.exp(1.077-1))**(-1)
    if z<1:
        c = 1
    else:
        c = ci + (1-ci1)
        
    SMF = HMF*np.gradient(np.log10(Mh_arr))/np.gradient(np.log10(Mstar_arr)) #varepsilon(Mh_arr,model_SFR,a)
    SMF = 10**(sigma_P(z)**2/2*np.log(10)*(np.gradient(np.log10(SMF))/np.gradient(np.log10(Mstar_arr)))**2)*SMF
    SMF = scipy.interpolate.interp1d(Mstar_arr,SMF, fill_value="extrapolate")
    SMF = f_passive_obs(Masses_star,a)*SMF(Masses_star*10**(-mu_SMF))+SMF(Masses_star*10**(-mu_SMF))*(1-f_passive_obs(Masses_star*10**(-kappa_SMF),a))
    SMF = c*SMF
    """
    return SMF
    

class NaNException(Exception):
    """Integrator hit a NaN."""

    pass


def hmf_integral_gtm(M, dndm, mass_density=False):
 
    # Eliminate NaN's
    m = M[np.logical_not(np.isnan(dndm))]
    dndm = dndm[np.logical_not(np.isnan(dndm))]
    dndlnm = m * dndm

    if len(m) < 4:
        raise NaNException(
            "There are too few real numbers in dndm: len(dndm) = %s, #NaN's = %s"
            % (len(M), len(M) - len(dndm))
        )

    # Calculate the mass function (and its integral) from the highest M up to 10**18
    if m[-1] < m[0] * 10**18 / m[3]:
        m_upper = np.arange(
            np.log(m[-1]), np.log(10**18), np.log(m[1]) - np.log(m[0])
        )
        mf_func = _spline(np.log(m), np.log(dndlnm), k=1)
        mf = mf_func(m_upper)

        if not mass_density:
            int_upper = intg.simps(np.exp(mf), dx=m_upper[2] - m_upper[1], even="first")
        else:
            int_upper = intg.simps(
                np.exp(m_upper + mf), dx=m_upper[2] - m_upper[1], even="first"
            )
    else:
        int_upper = 0

    # Calculate the cumulative integral (backwards) of [m*]dndlnm
    if not mass_density:
        ngtm = np.concatenate(
            (
                intg.cumtrapz(dndlnm[::-1], dx=np.log(m[1]) - np.log(m[0]))[::-1],
                np.zeros(1),
            )
        )
    else:
        ngtm = np.concatenate(
            (
                intg.cumtrapz(m[::-1] * dndlnm[::-1], dx=np.log(m[1]) - np.log(m[0]))[
                    ::-1
                ],
                np.zeros(1),
            )
        )

    return ngtm + int_upper

def Mh_EPS(z, k, Pk, rhoM, model_H, model, par1, par2, Mh0):
    zf = -0.0064*(np.log10(Mh0))**2+0.0237*np.log10(Mh0) + 1.8837
    q = 4.137*zf**(-0.9476)
    c_ST = 3.3
    R_M0 = (3.0*Mh0/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)
    R_M0q = (3.0*(Mh0/q)/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)

    func_EPS = 1/(np.sqrt(sigma(k,Pk,R_M0q)**2-sigma(k,Pk,R_M0)**2))
    deltac = delta_c_at_ac(model_H, model, 1, par1, par2)
    lineardelta_z0 = linear(1e-5, 1, model_H, model, par1, par2)[-1,1]
    lineardelta = scipy.interpolate.interp1d(linear(1e-5, 1, model_H, model, par1, par2)[:,0],linear(1e-5, 1, model_H, model, par1, par2)[:,1]/lineardelta_z0, fill_value = 'extrapolate')
    z_arr = np.linspace(0,1,50)
    dlineardz = np.gradient(lineardelta(1/(1+z_arr)))/np.gradient(z_arr)
    dlineardz_interp = scipy.interpolate.interp1d(z_arr, dlineardz, fill_value = 'extrapolate')
    dlineardz0 = dlineardz_interp(0)
    alpha = (deltac*np.sqrt(2/np.pi)*dlineardz0+1)*func_EPS
    beta = -func_EPS
    H = H_f(model_H,1/(1+z), par1, par2)
    Mh_EPS = Mh0*(1+z)**alpha*np.exp(beta*z)
    return [Mh_EPS, 71.6*(Mh0/1e12)*(h/0.7)*func_EPS*((1+z)-alpha/func_EPS)*H/(h*100)] #[EPS Mass, EPS Mass temporal derivative]

# SFR to MUV confertion, taken from https://github.com/XuejianShen/highz-empirical-variability
def convert_sfr_to_Muv(sfr, model_Muv):
    if model_Muv == 'Kennicutt2012':
        logCx = 43.35 # Kennicutt & Evans 2012 (assuming Kroupa IMF; using STARBURST99; solar metallicity)
        logLx = np.log10(sfr) + logCx  # log Lx in erg/s
        fnu = 10**logLx / (4*np.pi* (10*con.pc.value*100)**2 ) / (con.c.value/(1500*1e-10))
        Muv = -2.5 * np.log10(fnu) - 48.6 # AB mag
    elif model_Muv == "Madau2014":
        fnu = (sfr / 1.15e-28)/ (4*np.pi* (10*con.pc.value*100)**2 ) # erg/s/Hz/cm^2
        Muv = -2.5 * np.log10(fnu) - 48.6 # AB mag
    return Muv

# dust attenuation models, taken from https://github.com/XuejianShen/highz-empirical-variability 
def dust_attenuation(muv, dust_norm = "fixed"):
    # muv: intrinsic UV magnitude
    k_softplus = 10
    if dust_norm == "fixed":
        C0, C1 = 4.43, 1.99  # IRX beta relation, M99
        slope = -0.17; Mref = -19.5; intercept = -2.085 # Cullen 2023
        scatter = 0 # for a median relation
        #scatter=0.35 # for a mean relation

        prefactor = 1/(1 - C1 * slope)
        muv_obs = prefactor * (  muv  + C0 + C1 * intercept - C1 * slope * Mref + 0.2 * np.log(10) * C1**2 * scatter**2  )    # Vogelsberger 2020
        #return muv_obs * (muv_obs >= muv) + muv * (muv_obs < muv)
        return 1/k_softplus * np.log(1 + np.exp( k_softplus *( muv_obs - muv) )) + muv
    else:
        A = 1/(1 - 1.99 * (-0.17)) - 1 # Vogelsberger 2020
        B = dust_norm
        Auv = A * (muv - (-21)) + B
        muv_obs = muv + Auv
        #return muv_obs * (muv_obs >= muv) + muv * (muv_obs < muv)
        return 1/k_softplus * np.log(1 + np.exp( k_softplus *( muv_obs - muv) )) + muv    

