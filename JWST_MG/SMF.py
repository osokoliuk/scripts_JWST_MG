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

    # Code for Behroozi+13 SMHR is taken from
    # https://github.com/dr-guangtou/asap/tree/master
    def behroozi13_evolution(self, redshift):
        scale = 1.0 / (1.0 + redshift)
        scale_minus_one = -redshift / (1.0 + redshift)

        mh_1_0, mh_1_a, mh_1_z = 11.514, -1.793, -0.251
        epsilon_0, epsilon_a = -1.777, -0.006
        epsilon_z, epsilon_a_2 = -0.000, -0.119
        alpha_0, alpha_a = -1.412, 0.731
        delta_0, delta_a, delta_z = 3.508, 2.608, -0.043
        gamma_0, gamma_a, gamma_z = 0.316, 1.319, 0.279

        nu_a = np.exp(-4.0 * (scale ** 2.0))

        mh_1 = mh_1_0 + ((mh_1_a * scale_minus_one) + mh_1_z * redshift) * nu_a
        epsilon = epsilon_0 + ((epsilon_a * scale_minus_one) + epsilon_z * redshift) + epsilon_a_2 * scale_minus_one
        alpha = alpha_0 + (alpha_a * scale_minus_one) * nu_a
        delta = delta_0 + (delta_a * scale_minus_one + delta_z * redshift) * nu_a
        gamma = gamma_0 + (gamma_a * scale_minus_one + gamma_z * redshift) * nu_a

        return mh_1, epsilon, alpha, delta, gamma


    def behroozi13_f(self, x, alpha, delta, gamma):
        term_1 = -1.0 * np.log10(10.0 ** (alpha * x) + 1.0)

        term_2 = delta * (np.log10(1.0 + np.exp(x)) ** gamma) / (1.0 + np.exp(10.0 ** -x))

        return term_1 + term_2

    def behroozi13_mh_to_ms(self, halo_mass, z, **kwargs):
        logmh = np.log10(halo_mass)
        mh_1, epsilon, alpha, delta, gamma = self.behroozi13_evolution(z, **kwargs)

        mhalo_ratio = logmh - mh_1

        star_mass = mh_1 + epsilon + (self.behroozi13_f(mhalo_ratio, alpha, delta, gamma) -
                                self.behroozi13_f(0.0, alpha, delta, gamma))
        star_mass = 10**star_mass
        epsilon_star = star_mass/(halo_mass*Omegab0/Omegam0)
        return epsilon_star
        
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
            epstar = self.behroozi13_mh_to_ms(Mh, z)
        elif model_SFR == 'double_power':
            Mp = 10**12.1
            alo = -1.32
            ahi = 0.43
            f0 = 10**(-1.69)
            epstar = f0/((Mh/Mp)**alo + (Mh/Mp)**ahi)
            epstar = epstar/(Omegab0/Omegam0)
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

        z = 1/a-1
        mu_SMF = -0.020+0.081*(a-1)
        kappa_SMF = 0.045 + (-0.155)*(a-1)
        ci = 0.273*(1+np.exp(1.077-z))**(-1)
        ci1 = 0.273*(1+np.exp(1.077-1))**(-1)
        if z<1:
            c = 1
        else:
            c = ci + (1-ci1)
            
        
        #SMF = 10**(self.sigma_P(z)**2/2*np.log(10)*(np.gradient(np.log10(SMF))/np.gradient(np.log10(Masses_star)))**2)*SMF
        SMF = scipy.interpolate.interp1d(Masses_star,SMF, fill_value="extrapolate")
        SMF = self.f_passive_obs(Masses_star,a)*SMF(Masses_star*10**(-mu_SMF))+SMF(Masses_star*10**(-mu_SMF))*(1-self.f_passive_obs(Masses_star*10**(-kappa_SMF),a))
        SMF = c*SMF
        
        return Masses_star, SMF
