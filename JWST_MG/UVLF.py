from JWST_MG.constants import *
from JWST_MG.cosmological_functions import cosmological_functions
from JWST_MG.delta_c import delta_c
from JWST_MG.HMF import HMF
from JWST_MG.SMF import SMF
from JWST_MG.SMD import SMD

class UVLF:
   ########################################################################
    # Initialize a class UVLF (Ultra-Violet Luminosity Function)
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

    def Mh_EPS(self, a, rhoM, model_H, model, par1, par2, Mh0):
        Pk_library = HMF(a, model, model_H, par1, par2, self.Masses)
        cosmological_library = cosmological_functions(a, model, model_H, par1, par2)
        deltac_library = delta_c(a, model, model_H, par1, par2)
        
        z = 1/a-1
        c_ST = 3.3
        k = kvec/h
        zf = -0.0064*(np.log10(Mh0))**2+0.0237*np.log10(Mh0) + 1.8837
        q = 4.137*zf**(-0.9476)
        R_M0 = (3.0*Mh0/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)
        R_M0q = (3.0*(Mh0/q)/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)

        Pk = np.array(Pk_library.Pk(1,model,par1,par2))*h**3
        func_EPS = 1/(np.sqrt(Pk_library.sigma(k,Pk,R_M0q)**2-Pk_library.sigma(k,Pk,R_M0)**2))
        deltac = deltac_library.delta_c_at_ac(1, model, model_H, par1, par2)
        lineardelta_z0 = deltac_library.linear(1e-5, 1, model_H, model, par1, par2)[-1,1]
        lineardelta = scipy.interpolate.interp1d(deltac_library.linear(1e-5, 1, model_H, model, par1, par2)[:,0],deltac_library.linear(1e-5, 1, model_H, model, par1, par2)[:,1]/lineardelta_z0, fill_value = 'extrapolate')
        z_arr = np.linspace(0,1,50)
        dlineardz = np.gradient(lineardelta(1/(1+z_arr)))/np.gradient(z_arr)
        dlineardz_interp = scipy.interpolate.interp1d(z_arr, dlineardz, fill_value = 'extrapolate')
        dlineardz0 = dlineardz_interp(0)
        
        alpha = (deltac*np.sqrt(2/np.pi)*dlineardz0+1)*func_EPS
        beta = -func_EPS
        H = cosmological_library.H_f(a, model_H, par1, par2)
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

