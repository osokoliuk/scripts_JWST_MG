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

    def calculate_alphabeta(self, a, rhoM, model_H, model, par1, par2, Mass):
        Pk_library = HMF(a, model, model_H, par1, par2, self.Masses)
        cosmological_library = cosmological_functions(
            a, model, model_H, par1, par2)
        deltac_library = delta_c(a, model, model_H, par1, par2)

        z = 1/a-1
        c_ST = 3.3
        k = kvec/h
        deltac = deltac_library.delta_c_at_ac(1, model, model_H, par1, par2)
        lineardelta_z0 = deltac_library.linear(
            1e-5, 1, model_H, model, par1, par2)[-1, 1]
        lineardelta = scipy.interpolate.interp1d(deltac_library.linear(1e-5, 1, model_H, model, par1, par2)[
                                                 :, 0], deltac_library.linear(1e-5, 1, model_H, model, par1, par2)[:, 1]/lineardelta_z0, fill_value='extrapolate')
        z_arr = np.linspace(0, 1, 50)
        dlineardz = np.gradient(lineardelta(1/(1+z_arr)))/np.gradient(z_arr)
        dlineardz_interp = scipy.interpolate.interp1d(
            z_arr, dlineardz, fill_value='extrapolate')
        dlineardz0 = dlineardz_interp(0)
        alpha = []
        beta = []
        Masses = np.logspace(8, 16, 500)
        Pk = np.array(Pk_library.Pk(1, model, par1, par2))*h**3

        for Mh0 in Masses:
            zf = -0.0064*(np.log10(Mh0))**2+0.0237*np.log10(Mh0) + 1.8837
            q = 4.137*zf**(-0.9476)
            R_M0 = (3.0*Mh0/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)
            R_M0q = (3.0*(Mh0/q)/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)

            func_EPS = 1/(np.sqrt(Pk_library.sigma(k, Pk, R_M0q)
                          ** 2-Pk_library.sigma(k, Pk, R_M0)**2))

            alpha.append((deltac*np.sqrt(2/np.pi)*dlineardz0+1)*func_EPS)
            beta.append(-func_EPS)

        alpha = scipy.interpolate.interp1d(
            Masses, alpha, fill_value='extrapolate')
        beta = scipy.interpolate.interp1d(
            Masses, beta, fill_value='extrapolate')
        return alpha(Mass), beta(Mass)

    def MAR(self, a, rhoM, model_H, model, par1, par2, Masses):
        z = 1/a-1
        Pk_library = HMF(a, model, model_H, par1, par2, Masses)
        cosmological_library = cosmological_functions(
            a, model, model_H, par1, par2)
        deltac_library = delta_c(a, model, model_H, par1, par2)

        alpha, beta = self.calculate_alphabeta(
            a, rhoM, model_H, model, par1, par2, Masses)

        if hasattr(a, '__len__') and (not isinstance(a, str)):
            dMh_EPS = []
            for ai in a:
                zi = 1/ai-1
                H = cosmological_library.H_f(ai, model_H, par1, par2)
                Mh_EPS = Masses*(1+zi)**alpha*np.exp(beta*zi)
                dMh_EPS.append(71.6*(Masses/1e12)*(h/0.7) *
                               (-alpha-beta*(1+zi))*H/(h*100))
        else:
            H = cosmological_library.H_f(a, model_H, par1, par2)
            Mh_EPS = Masses*(1+z)**alpha*np.exp(beta*z)
            dMh_EPS = 71.6*(Masses/1e12)*(h/0.7)*(-alpha-beta*(1+z))*H/(h*100)

        return dMh_EPS

    def SFR(self, a, rhoM, model, model_H, model_SFR, par1, par2, Masses, f0):
        dMdt = self.MAR(a, rhoM, model_H, model, par1, par2, Masses)
        SMF_library = SMF(a, model, model_H, model_SFR, par1, par2, Masses, f0)
        fstar = SMF_library.epsilon(Masses, model_SFR, a, f0)*Omegab0/Omegam0
        SFR = fstar*dMdt
        return SFR

    def SFRD(self, a, rhoM, model, model_H, model_SFR, par1, par2, Masses, f0):
        SFR_fid = self.SFR(a, rhoM, model, model_H,
                           model_SFR, par1, par2, Masses, f0)
        SMD_library = SMD(a, model, model_H, model_SFR, par1, par2, Masses, f0)
        HMF_library = HMF(a, model, model_H, par1, par2, Masses)
        HMF_fid = HMF_library.ST_mass_function(
            rhoM, Masses, a, model_H, model, par1, par2)

        SFRD = SMD_library.sfr_integral_gtm(
            Masses, HMF_fid, SFR_fid, mass_density=True)
        return SFRD

    # SFR to MUV confertion
    # Taken from https://github.com/XuejianShen/highz-empirical-variability
    def convert_sfr_to_Muv(self, sfr, model_Muv="Madau2014"):
        if model_Muv == 'Kennicutt2012':
            # Kennicutt & Evans 2012 (assuming Kroupa IMF; using STARBURST99; solar metallicity)
            logCx = 43.35
            logLx = np.log10(sfr) + logCx  # log Lx in erg/s
            fnu = 10**logLx / (4*np.pi * (10*con.pc.value*100)
                               ** 2) / (con.c.value/(1500*1e-10))
            Muv = -2.5 * np.log10(fnu) - 48.6  # AB mag
        elif model_Muv == "Madau2014":
            fnu = (sfr / 1.15e-28) / (4*np.pi *
                                      (10*con.pc.value*100)**2)  # erg/s/Hz/cm^2
            Muv = -2.5 * np.log10(fnu) - 48.6  # AB mag
        return Muv

    # Account for dust attenuation in MUV
    # Taken from https://github.com/XuejianShen/highz-empirical-variability
    def dust_attenuation(self, muv, dust_norm="fixed"):
        # muv: intrinsic UV magnitude
        k_softplus = 10
        if dust_norm == "fixed":
            C0, C1 = 4.43, 1.99  # IRX beta relation, M99
            slope = -0.17
            Mref = -19.5
            intercept = -2.085  # Cullen 2023
            scatter = 0  # for a median relation
            # scatter=0.35 # for a mean relation

            prefactor = 1/(1 - C1 * slope)
            muv_obs = prefactor * (muv + C0 + C1 * intercept - C1 * slope *
                                   Mref + 0.2 * np.log(10) * C1**2 * scatter**2)    # Vogelsberger 2020
            # return muv_obs * (muv_obs >= muv) + muv * (muv_obs < muv)
            return 1/k_softplus * np.log(1 + np.exp(k_softplus * (muv_obs - muv))) + muv
        else:
            A = 1/(1 - 1.99 * (-0.17)) - 1  # Vogelsberger 2020
            B = dust_norm
            Auv = A * (muv - (-21)) + B
            muv_obs = muv + Auv
            # return muv_obs * (muv_obs >= muv) + muv * (muv_obs < muv)
            return 1/k_softplus * np.log(1 + np.exp(k_softplus * (muv_obs - muv))) + muv

    # Map SFR to MUV with dust correction
    # Taken from https://github.com/XuejianShen/highz-empirical-variability
    def mapfunc_mhalo_to_muv(self, a, rhoM, model, model_H, model_SFR, par1, par2, Masses, f0, dust_norm="fixed", include_dust=True):
        '''
        mapping from halo mass to UV magnitude (without scatter)
        muv: UV magnitude
        log_mhalo: log10 of halo mass in Msun
        '''
        sfr = self.SFR(a, rhoM, model, model_H,
                       model_SFR, par1, par2, Masses, f0)

        muv_raw = self.convert_sfr_to_Muv(sfr)
        if include_dust:
            muv = self.dust_attenuation(muv_raw, dust_norm=dust_norm)
        else:
            muv = muv_raw
        muv = scipy.interpolate.interp1d(Masses, muv, fill_value='extrapolate')
        return muv(Masses)

    # Derive factor dMUV/dlogMh
    # Taken from https://github.com/XuejianShen/highz-empirical-variability

    def mapfunc_jacobian_numeric(self, a, rhoM, model, model_H, model_SFR, par1, par2, Masses, f0, dust_norm="fixed", include_dust=True):
        muv = self.mapfunc_mhalo_to_muv(
            a, rhoM, model, model_H, model_SFR, par1, par2, Masses, f0, dust_norm="fixed", include_dust=True)
        dmuv_dlogm = np.gradient(muv)/np.gradient(Masses)
        return np.abs(dmuv_dlogm)

    # Convolve MUV with Gaussian kernel of width sigma_uv
    # This takes into account scatter that may arise from different factors,
    # For example from the SMHR uncertainties and HMF choice
    # Taken from https://github.com/XuejianShen/highz-empirical-variability
    def convolve_on_grid(self, input_grid, input_weight, sigma_uv):
        grid_binsize = input_grid[1] - input_grid[0]
        minimum_sigma = grid_binsize/4.  # set to the binsize divided by a constant
        # regulate the miminum sigma to be of order the binsize (~ 0.01 dex)
        sigma_uv = max(sigma_uv, minimum_sigma)

        output_weight = np.zeros(len(input_grid))
        for i, mapfrom in enumerate(input_grid):
            raw_output = np.zeros(len(input_grid))
            for j, mapto in enumerate(input_grid):
                raw_output[j] += 1./np.sqrt(2*np.pi*sigma_uv**2) * \
                    np.exp(-0.5 * (mapto - mapfrom)**2 / sigma_uv**2)
            sum_raw_output = np.sum(raw_output)
            if sum_raw_output > 0:
                raw_output = raw_output/sum_raw_output
            else:
                raw_output = np.zeros(len(input_grid))
            output_weight += input_weight[i] * raw_output
        return output_weight

    # Finally compute UVLF by using all of the previously defined functions in this class
    # Taken from https://github.com/XuejianShen/highz-empirical-variability
    def compute_uv_luminosity_function(self, a, rhoM, model, model_H, model_SFR, par1, par2, Masses, f0, sigma_uv, dust_norm="fixed", include_dust=True):
        HMF_library = HMF(a, model, model_H, par1, par2, Masses)
        phi_halo_arr = HMF_library.ST_mass_function(
            rhoM, Masses, a, model_H, model, par1, par2)
        muv_arr = self.mapfunc_mhalo_to_muv(
            a, rhoM, model, model_H, model_SFR, par1, par2, Masses, f0, dust_norm, include_dust)
        dmuv_dlogm = self.mapfunc_jacobian_numeric(
            a, rhoM, model, model_H, model_SFR, par1, par2, Masses, f0, dust_norm, include_dust)
        if sigma_uv > 0:
            phi_uv_arr = self.convolve_on_grid(
                muv_arr, phi_halo_arr/dmuv_dlogm, sigma_uv=sigma_uv)
        else:
            phi_uv_arr = phi_halo_arr/dmuv_dlogm
        # print(phi_uv_arr)
        return muv_arr, phi_uv_arr
