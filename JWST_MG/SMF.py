from JWST_MG.constants import *
from JWST_MG.cosmological_functions import cosmological_functions
from JWST_MG.delta_c import delta_c
from JWST_MG.HMF import HMF

from hmf import MassFunction     # The main hmf class
import numpy
from numpy import log, log10, sqrt, exp, pi, arange, logspace, linspace, r_, c_
# Don't judge me: 
from matplotlib.pyplot import *
from scipy.interpolate import UnivariateSpline as US
import hmf
from astropy.cosmology import FlatLambdaCDM

class SMF:
   ########################################################################
    # Initialize a class SMF (Stellar Mass Function)
    # float a - scale factor value (related to redshift via a = 1/(1+z))
    # string model - model of MG for the derivation of mu parameter
    # string model_H - model of MG for H(a)
    # string model_SFR - model that determines SMHR
    # float par1, par2 - corresponding MG parameters
    # Masses - array of CDM halo masses
    # f0 - parameter for "double-power" SMHR, otherwise None
    ########################################################################

    def __init__(self, a, model, model_H, model_SFR, par1, par2, Masses, f0=None):
        self.a = a
        self.model = model
        self.model_H = model_H
        self.model_SFR = model_SFR
        self.par1 = par1
        self.par2 = par2
        self.Masses = Masses
        self.f0 = f0

    def func_SFR(self, x, a):
        z = 1/a-1
        nu = np.exp(-4*a**2)
        alpha_SFR = -1.412+0.731*(a-1)*nu
        Delta_SFR = 3.508 + (2.608*(a-1)-0.043*z)*nu
        gamma_SFR = 0.316+(1.319*(a-1)+0.279*z)*nu
        return -np.log10(10**(alpha_SFR*x)+1)+Delta_SFR*(np.log10(1+np.exp(x)))**gamma_SFR/(1+np.exp(10**(-x)))

    def epsilon(self, Mh, model_SFR, a, f0):
        z = 1/a-1
        if model_SFR == 'toy':
            epstar = f0
        elif model_SFR == 'phenomenological_regular':
            if z < 10:
                epstar = 0.15 - 0.03*(z-6)
            else:
                epstar = 0.03
        elif model_SFR == 'phenomenological_extreme':
            epstar = 1
        elif model_SFR == 'Behroozi':
            nu = np.exp(-4*a**2)
            log10M1 = 11.514+(-1.793*(a-1)+(-0.251)*z)*nu
            log10eps = -1.777+(-0.006*(a-1)+(-0.000)*z)*nu-0.119*(a-1)
            log10Mstar = log10eps + log10M1 + \
                self.func_SFR(np.log10(Mh)-log10M1, a)-self.func_SFR(0, a)
            Mstar = 10**log10Mstar
            epstar = (Mstar/Mh)/(Omegab0/Omegam0)
        elif model_SFR == 'double_power':
            Mp = 2.8*10**11
            alo = 0.49
            ahi = -0.61
            epstar = 2*f0/((Mh/Mp)**alo + (Mh/Mp)**ahi)
        else:
            raise Exception("Incorrect SFR model used.")
        return epstar

    def sigma_P(self, z):
        sigma0 = 0.07
        sigmaz = 0.05
        return sigma0 + sigmaz*z

    def f_passive_obs(self, Masses_star, a):
        z = 1/a-1
        return ((Masses_star/(10**(10.2+0.5*z)))**(-1.3)+1)**(-1)

    def SMF_obs(self, Masses, rhoM, a, model_H, model, model_SFR, par1, par2, k, Pk, f0):
        # Mh_arr = np.logspace(6.5,18,1000)
        HMF_library = HMF(a, model, model_H, par1, par2, Masses)
        HMF_fid = HMF_library.ST_mass_function(
            rhoM, Masses, a, model_H, model, par1, par2, k, Pk)
        #HMF_fid = np.log(10)*Masses*HMF_fid

        Masses_star = self.epsilon(
            Masses, model_SFR, a, f0)*Omegab0/Omegam0*Masses
        # varepsilon(Mh_arr,model_SFR,a)
        SMF = Masses_star*np.log(10)*HMF_fid*np.gradient(Masses)/ \
            np.gradient(Masses_star)

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
        return Masses_star, SMF
