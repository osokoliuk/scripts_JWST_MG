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

    # Rodriguez-Puebla parameterization is taken from
    # the repo https://github.com/dr-guangtou/asap/blob/master/asap/shmr.py
    def puebla17_p(self, x, y, z):
        """The P(x, y, z) function used in Rodriguez-Puebla+17."""
        return y * z - (x * z) / (1.0 + z)

    def puebla17_q(self, z):
        """The Q(z) function used in Rodriguez-Puebla+17."""
        return np.exp(-4.0 / (1.0 + z) ** 2.0)

    def puebla17_g(self, x, alpha, delta, gamma):
        """The g(x) function used in Behroozi+13."""
        term_1 = -np.log10(10.0 ** (-alpha * x) + 1.0)

        term_2 = delta * (np.log10(1.0 + np.exp(x)) ** gamma) / (1.0 + np.exp(10.0 ** -x))

        return term_1 + term_2

    def puebla17_evolution(self, redshift):
        """Parameterize the evolution in term of scale factor.

        Using the best-fit parameters in Rodriguez-Puebla+17.
        """

        mh_1_0, mh_1_1, mh_1_2 = 11.548, -1.297, -0.026
        epsilon_0, epsilon_1 = -1.758, 0.110
        epsilon_2, epsilon_3 = -0.061, -0.023
        alpha_0, alpha_1, alpha_2 = 1.975, 0.714, 0.042
        delta_0, delta_1, delta_2 = 3.390, -0.472, -0.931
        gamma_0, gamma_1 = 0.498, -0.157

        mh_1 = mh_1_0 + self.puebla17_p(mh_1_1, mh_1_2, redshift) * self.puebla17_q(redshift)
        epsilon = epsilon_0 + (self.puebla17_p(epsilon_1, epsilon_2, redshift) * self.puebla17_q(redshift) +
                            self.puebla17_p(epsilon_3, 0.0, redshift))
        alpha = alpha_0 + self.puebla17_p(alpha_1, alpha_2, redshift) * self.puebla17_q(redshift)
        delta = delta_0 + self.puebla17_p(delta_1, delta_2, redshift) * self.puebla17_q(redshift)
        gamma = gamma_0 + self.puebla17_p(gamma_1, 0.0, redshift) * self.puebla17_q(redshift)

        return mh_1, epsilon, alpha, delta, gamma

    def puebla17_mh_to_ms(self, Mh, z):
        logmh = np.log10(Mh)
        mh_1, epsilon, alpha, delta, gamma = self.puebla17_evolution(z)

        mhalo_ratio = logmh - mh_1

        Mstar = mh_1 + epsilon + (self.puebla17_g(mhalo_ratio, alpha, delta, gamma) -
                                self.puebla17_g(0.0, alpha, delta, gamma))
        Mstar = 10**Mstar
        epsilon_star = Mstar/(Mh*Omegab0/Omegam0)
        return epsilon_star
        
    def epsilon(self, Mh, model_SFR, a, f0):
        z = 1/a-1
        if model_SFR == 'phenomenological_extreme':
            epstar = 1
        elif model_SFR == 'Puebla':
            epstar = self.puebla17_mh_to_ms(Mh, z)
        elif model_SFR == 'double_power':
            Mp = 10**11.16
            alo = 0.8
            ahi = -0.53
            f0 = 10**(-1.26)
            C10 = (1e10/Mp)**(-alo) + (1e10/Mp)**(-ahi)
            epstar = C10*f0/((Mh/Mp)**(-alo) + (Mh/Mp)**(-ahi))
        else:
            raise Exception("Incorrect SFR model used.")
        return epstar

    def sigma_P(self, z):
        sigma0 = 0.1
        sigmaz = 0.05
        return sigma0 + sigmaz*z

    def G_prob(self, x, Mstar, z):
        sigma = self.sigma_P(z)
        return 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-1/(2*sigma**2)*np.log10(Mstar/x)**2)

    def f_passive_obs(self, Masses_star, a):
        z = 1/a-1
        return ((Masses_star/(10**(10.2+0.5*z)))**(-1.3)+1)**(-1)

    def SMF_interpolation(self, Mstar_var, SMF, z, down, up):
        return scipy.integrate.quad(lambda logx: self.G_prob(10**logx, Mstar_var, z)*SMF(10**logx), down, up)[0]

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

        if model_SFR == 'Puebla':
            z = 1/a-1
            mu_SMF = -0.020+0.081*(a-1)
            kappa_SMF = 0.045 + (-0.155)*(a-1)
            ci = 0.273*(1+np.exp(1.077-z))**(-1)
            ci1 = 0.273*(1+np.exp(1.077-1))**(-1)
            if z<1:
                c = 1
            else:
                c = ci + (1-ci1)

            SMF = scipy.interpolate.interp1d(Masses_star,SMF, fill_value="extrapolate")
            Mstar_grid = np.logspace(6,12.1,50)
            SMF_arr = []
            for Mstar in Mstar_grid:
                SMF_arr.append(self.SMF_interpolation(Mstar, SMF, z,np.log10(min(Masses_star)),np.log10(max(Masses_star))))
            SMF = c*np.array(SMF_arr)
            Masses_star = Mstar_grid

            #Masses_star = 10**mu_SMF*Masses_star 

        """
        SMF = 10**(self.sigma_P(z)**2/2*np.log(10)*(np.gradient(np.log10(SMF))/np.gradient(np.log10(Masses_star)))**2)*SMF
        SMF = scipy.interpolate.interp1d(Masses_star,SMF, fill_value="extrapolate")
        SMF = self.f_passive_obs(Masses_star,a)*SMF(Masses_star*10**(-mu_SMF))+SMF(Masses_star*10**(-mu_SMF))*(1-self.f_passive_obs(Masses_star*10**(-kappa_SMF),a))
        SMF = c*SMF
        """
        return Masses_star, SMF
