from JWST_MG.constants import *
from JWST_MG.cosmological_functions import cosmological_functions
from JWST_MG.delta_c import delta_c


class HMF:
    ########################################################################
    # Initialize a class HMF (Halo Mass Function)
    # float a - scale factor value (related to redshift via a = 1/(1+z))
    # array of floats k - wavenumber in units of 1/Mpc
    # string model - model of MG for the derivation of mu parameter
    # string model_H - model of MG for H(a)
    # float par1, par2 - corresponding MG parameters
    # Masses - array of CDM halo masses
    ########################################################################

    def __init__(self, a, model, model_H, par1, par2, Masses):
        self.a = a
        self.model = model
        self.model_H = model_H
        self.par1 = par1
        self.par2 = par2
        self.Masses = Masses
        self.common_settings = {'n_s': 0.9665,
                           'A_s': 2.101e-9,
                           'tau_reio': 0.0561,
                           'omega_b': 0.02242,
                           'omega_cdm': 0.11933,
                           'h': 0.6766,
                           'YHe': 0.2425,
                           'T_cmb': 2.7255,
                           'gauge': 'newtonian',  # FOR MGCLASS TO WORK, GAUGE NEEDS TO BE NEWTONIAN
                           'k_pivot': 0.05,
                           'mg_z_init': 99.000,
                           'P_k_max_1/Mpc': 10.0,
                           'perturb_sampling_stepsize': 0.05,
                           'output': 'mPk',
                           'z_max_pk': 99}
        self.M = Class()

    def Pk(self, a, model, par1, par2):
        common_settings = self.common_settings
        common_settings['mg_ansatz'] = model
        if model == 'plk_late':
            common_settings['mg_E11'] = par1
            common_settings['mg_E22'] = par2
        elif model == 'z_flex_late':
            common_settings['mg_muz'] = par1
            common_settings['mg_gamz'] = par2
            common_settings['mg_zzn'] = 1
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
        elif model == 'kmoufl':
            common_settings['beta_kmfl'] = par1
            common_settings['k0_kmfl'] = par2
        else:
            raise Exception("The chosen model is not recognised")
        M = self.M
        M.set(common_settings)
        M.compute()

        Pk = []
        for k in kvec:
            Pk.append(M.pk(k, 1/a-1))
        M.empty()
        M.struct_cleanup()
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

    def sigma(self, k, Pk, R):
        yinit = np.array([0.0], dtype=np.float64)
        eps = 1e-13  # change this for higher/lower accuracy
        h1 = 1e-12
        hmin = 0.0
        beta_ST = 4.8
        W = (1+(k*R)**beta_ST)**(-1)
        Pk1 = Pk*W**2*k**2/(2.0*np.pi**2)

        return np.sqrt(IL.odeint(yinit, k[0], k[-1], eps,
                                 h1, hmin, np.log10(k), Pk1,
                                 'sigma', verbose=False)[0])

    def dSdM(self, k, Pk, rhoM, M):
        c_ST = 3.3
        R1 = (3.0*M/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)
        s1 = self.sigma(k, Pk, R1)

        M2 = M*1.0001
        R2 = (3.0*M2/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)
        s2 = self.sigma(k, Pk, R2)

        return (s2-s1)/(M2-M)

    # Taken from https://pylians3.readthedocs.io/en/master/mass_function.html
    # And properly modified to incorporate MG theories with varying delta_c
    def ST_mass_function(self, rhoM, Masses, a, model_H, model, par1, par2):
        c_ST = 3.3
        deltac = delta_c(a, model, model_H, par1, par2)
        deltac = deltac.delta_c_at_ac(a, model, model_H, par1, par2)
        Pk = np.array(self.Pk(a, model, par1, par2))*h**3
        k = kvec/h
        if hasattr(Masses, '__len__') and (not isinstance(Masses, str)):
            dndM = np.zeros(Masses.shape[0], dtype=np.float64)
            for i, M in enumerate(Masses):
                R = (3.0*M/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)
                nu = (deltac/self.sigma(k, Pk, R))**2
                dndM[i] = -(rhoM/M)*self.dSdM(k, Pk, rhoM, M) / \
                    self.sigma(k, Pk, R)
                dndM[i] *= 0.3222*np.sqrt(2*nu/np.pi)*(1+1/(nu**0.3))
                dndM[i] *= np.exp(-0.5*nu)
        else:
            R = (3.0*Masses/(4.0*np.pi*rhoM*c_ST**3))**(1.0/3.0)
            nu = (deltac/self.sigma(k, Pk, R))**2

            dndM = -(rhoM/Masses)*self.dSdM(k, Pk,
                                            rhoM, Masses)/self.sigma(k, Pk, R)
            dndM *= 0.3222*np.sqrt(2*nu/np.pi)*(1+1/(nu**0.3))
            dndM *= np.exp(-0.5*nu)

        return dndM
