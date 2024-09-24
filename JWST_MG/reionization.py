from JWST_MG.constants import *
from JWST_MG.cosmological_functions import cosmological_functions
from JWST_MG.delta_c import delta_c
from JWST_MG.HMF import HMF
from JWST_MG.SMF import SMF
from JWST_MG.SMD import SMD
from JWST_MG.UVLF import UVLF

class reionization:
    ########################################################################
    # Initialize a class delta_c (critical threshold for linear density perturbations)
    # float ac - scale factor, at which spherical collapse happens
    # string model - model of MG for the derivation of mu parameter
    # string model_H - model of MG for H(a)
    # float par1, par2 - corresponding MG parameters
    # float delta_i - initial linear/non-linear overdensity at a = 1e-5
    ########################################################################

    def __init__(self, a_arr, model, model_H, par1, par2, model_SFR=None, f0=None):
        self.a_arr = a_arr
        self.ac = a_arr[-1]
        self.model = model
        self.model_H = model_H
        self.par1 = par1
        self.par2 = par2

        global cosmological_library
        cosmological_library = cosmological_functions(
            a_arr, model, model_H, par1, par2)
        global deltac_library
        deltac_library = delta_c(
            a_arr, model, model_H, par1, par2)
        self.deltai = deltac_library.binary_search_di(
            self.ac, self.model, self.model_H, self.par1, self.par2, 0, len(delta_ini), abs_err)
        self.delta_nl = deltac_library.non_linear(
            self.deltai, self.a_arr, self.model, self.model_H, self.par1, self.par2)
	
    def delta_nl_a(self, x):
        func = scipy.interpolate.interp1d(
            self.a_arr, self.delta_nl, fill_value="extrapolate")

        return func(x)


    def radius_evolution(self, y, a, model, model_H, par1, par2, a_arr):
        R, dRda = y
        cosmological_library = cosmological_functions(
        ai, model, model_H, par1, par2)
        delta_nl = self.delta_nl_a(a)
        H = cosmological_library.H_f(a, model_H, par1, par2)
        dH = cosmological_library.dH_f(a, model_H, par1, par2)
        Hprime = a*dH
        Rprime = a*dRda

        if model_H == "LCDM" or model_H == "wCDM":
            mu = cosmological_library.mu(a, model, model_H, par1, par2)
            ddRda = (-Hprime/H*Rprime +
                 (1+Hprime/H)*R - Omegam0*a**(-3)*H0**2 /
                 (2*H**2) * mu*(R+a/ai)*delta_nl - a*dRda)/a**2
        elif model_H == "nDGP":
            H_dot = a*H*dH
            beta = 1 + 2*H*par1/c*(1+H_dot/(3*H**2))
            epsilon = 8/(9*beta**2)*(H0*par1/c)**2*Omegam0*a**(-3)
            RRV = (epsilon*delta_nl)**(-1/3)
            mu = cosmological_library.mu(
                a, model, model_H, par1, par2, type='nonlinear', x=RRV)
            ddRda = (-Hprime/H*Rprime +
                     (1+Hprime/H)*R - Omegam0*a**(-3)*H0**2 /
                     (2*H**2) * mu*(R+a/ai)*delta_nl - a*dRda)/a**2
        elif model_H == "kmoufl":
            A_kmfl = 1.0 + par1*a
            X_kmfl = 0.5 * A_kmfl**2*(H*a)**2/((1-Omegam0-Omegar0)*H0**2)
            k_prime_mfl = 1.0 + 2.0*par2*X_kmfl
            epsl1_kmfl = 2.0*par1**2/k_prime_mfl
            epsl2_kmfl = a*par1/(1.0+par1*a)
            Hdot = a*H*dH
            ddRda = (H0**2*(1 + epsl1_kmfl)*Omegam0*(-1 + R**3) - 2*a**4*(Hdot + H**2*(3 + epsl2_kmfl))*R**2*dRda)/ (2.*a**5*H**2*R**2)
        return [dRda, ddRda]

    def radius_solve(self, model, model_H, par1, par2, a_arr):
        deltai = self.deltai
        cosmological_library = cosmological_functions(
            ai, model, model_H, par1, par2)
        Hi = cosmological_library.H_f(ai, model_H, par1, par2)
        if model != "kmoufl":
            R_arr = scipy.integrate.odeint(self.radius_evolution, [0, -ai*Hi*deltai/(3*(1+deltai))], a_arr, args=(
                model, model_H, par1, par2, a_arr), tfirst=False)[:, 0]

            return R_arr + a_arr/ai
        else:
            delta_nl = self.delta_nl_a(a_arr)
            R_arr = (a_arr**3*ai**6*(1+delta_nl)**2*(1+deltai))**(1/3)/(ai**3*(1+delta_nl))
            return R_arr
            #R_arr = scipy.integrate.odeint(self.radius_evolution, [(1+deltai)**(-1/3), -deltai/(3*ai*(1+deltai)**(4/3))], a_arr, args=(
            #    model, model_H, par1, par2, a_arr), tfirst=False)[:, 0]
            #return a_arr*R_arr

    def virial_theorem(self, model, model_H, par1, par2, a_arr):
        G = 1/(8*np.pi)
        ac = a_arr[-1]
        cosmological_library = cosmological_functions(
            ac, model, model_H, par1, par2)
        R_arr = self.radius_solve(
            model, model_H, par1, par2, a_arr)
        H_arr = cosmological_library.H_f(a_arr, model_H, par1, par2)
        dH_arr = cosmological_library.dH_f(a_arr, model_H, par1, par2)
        delta_nl = self.delta_nl_a(a_arr)

        R_arr[R_arr == -inf] = 0
        R_arr[R_arr == inf] = 0
        R_arr[R_arr < 0] = 0
        R_arr[a_arr < ai] = 0
        R_arr[a_arr > ac] = 0

        H_dot = a_arr*H_arr*dH_arr
        Rdot = a_arr*H_arr*np.gradient(R_arr)/np.gradient(a_arr)
        rho = 3*H0**2*Omegam0
        M = 4*np.pi/3*R_arr**3*rho*a_arr**(-3)*(1+delta_nl)
        deltaM = delta_nl/(1+delta_nl)*M
        T = 3/10*M*Rdot**2

        if model_H == "LCDM" or model_H == "wCDM" or model_H == "kmoufl":
            mu_arr = cosmological_library.mu(a_arr, model, model_H, par1, par2)
            U = -3/5*G*mu_arr*M*deltaM/R_arr+3/5*(H_dot+H_arr**2)*M*R_arr**2
            virial = T+1/2*U
        elif model_H == "nDGP":
            beta = 1 + 2*H_arr*par1/c*(1+H_dot/(3*H_arr**2))
            epsilon = 8/(9*beta**2)*(H0*par1/c)**2*Omegam0*a_arr**(-3)
            RRV = (epsilon*delta_nl)**(-1/3)
            mu_arr = cosmological_library.mu(
                a_arr, model, model_H, par1, par2, type='nonlinear', x=RRV)
            Geff = G*mu_arr
            DeltaGeff = Geff - G
            U = -3/5*G*M**2/R_arr - 3/5*DeltaGeff*M * \
                deltaM/R_arr+3/5*(H_dot+H_arr**2)*M*R_arr**2
            virial = T+1/2*U
        else:
            raise ValueError("Unknown cosmology type!")

        virial[virial < 0] = 1e60

        virial = scipy.interpolate.interp1d(
            a_arr, virial, fill_value='extrapolate')

        a_turn = a_arr[R_arr.argmax()]
        a_vir = scipy.optimize.minimize_scalar(
            lambda x: virial(x), bounds=(ai, self.ac), method="bounded").x

        if model != "nDGP":
            return a_turn, a_vir
        else:
            return a_turn, a_vir, a_arr, mu_arr

    def Delta_vir(self, model, model_H, par1, par2, a_arr):
        ac = a_arr[-1]

        if model != "nDGP":
            a_turn, a_vir = self.virial_theorem(
                model, model_H, par1, par2, a_arr)
            Deltavir = (1+self.delta_nl_a(a_vir))*(ac/a_vir)**3
            return a_vir, Deltavir

        else:
            a_turn, a_vir, a_arr, mu_arr = self.virial_theorem(
                model, model_H, par1, par2, a_arr)
            Deltavir = (1+self.delta_nl_a(a_vir))*(ac/a_vir)**3
            return a_vir, Deltavir, a_arr, mu_arr

    def minimum_Mhalo(self, model, model_H, par1, par2, a_arr):
        ac = a_arr[-1]
        cosmological_library = cosmological_functions(
            ac, model, model_H, par1, par2)
        if model != "nDGP":
            a_vir, Deltavir = self.Delta_vir(
                model, model_H, par1, par2, a_arr)
            mu = cosmological_library.mu(ac, model, model_H, par1, par2)
        else:
            a_vir, Deltavir, ac_arr, mu_arr  = self.Delta_vir(
                model, model_H, par1, par2, a_arr)
            mu = mu_arr[np.argmin(np.abs(ac_arr-ac))]

        H = cosmological_library.H_f(ac, model_H, par1, par2)
        dH = cosmological_library.dH_f(ac, model_H, par1, par2)
        H_dot = ac*H*dH

        Tmin = 10000
        Ceff = -3/(4*np.pi*GN*rho*ac**(-3))*(H_dot+H**2)
        Mhalo_min = 1e-9*(10*kB*Tmin/(3*mu_mol*mP))**(3/2)*(GN*H0*np.sqrt(Omegam0 *
                                                                     ac**(-3)*Deltavir/2))**(-1)*(1/Deltavir*(Ceff)+mu*(1-1/Deltavir))**(-3/2)

        return a_vir, Mhalo_min

    def n_ion(self, a, rhoM, model, model_H, model_SFR, par1, par2, k, Pk, f0=None):
        #Nion = 10**53.14
        a_arr = np.linspace(ai, a, 1000)
        a_vir, Mhalo_min = self.minimum_Mhalo(
            model, model_H, par1, par2, a_arr)
        Masses = np.logspace(np.log10(Mhalo_min), 18, 1000)

        UVLF_library = UVLF(a, model, model_H, model_SFR,
                            par1, par2, Masses, f0)
        SFRD_fid = UVLF_library.SFRD(
            a, rhoM, model, model_H, model_SFR, par1, par2, Masses, k, Pk, f0)
        
        nion = SFRD_fid

        return nion

    # This takes into account scatter that may arise from different factors,
    # For example from the SMHR uncertainties and HMF choice
    # Taken from https://github.com/XuejianShen/highz-empirical-variability

    def convolve_on_grid(self, input_grid, input_weight, sigma_uv_input, normalization_check=False, verbose=False):
        # convolve with a gaussian kernel
        # over the UV luminosity function
        grid_binsize  = np.abs(input_grid[2:] - input_grid[:-2])/2.
        grid_binsize  = np.append(grid_binsize, grid_binsize[-1])
        grid_binsize  = np.append(grid_binsize[0], grid_binsize)

        minimum_sigma = 0.2 # set to the binsize divided by a constant
        #print(minimum_sigma)
        if np.isscalar(sigma_uv_input):
            sigma_uv = max(sigma_uv_input, minimum_sigma) # regulate the miminum sigma to be of order the binsize (~ 0.01 dex)
        else:
            sigma_uv = sigma_uv_input.copy()
            sigma_uv[sigma_uv < minimum_sigma] = minimum_sigma

        output_weight = np.zeros(len(input_grid))
        for i, mapfrom in enumerate(input_grid):
            raw_output    = np.zeros(len(input_grid))
            for j, mapto in enumerate(input_grid):
                if np.isscalar(sigma_uv):
                    raw_output[j] += 1./np.sqrt(2*np.pi*sigma_uv**2) * np.exp( -0.5 * (mapto - mapfrom)**2 / sigma_uv**2 ) * grid_binsize[j]
                else:
                    raw_output[j] += 1./np.sqrt(2*np.pi*sigma_uv[i]**2) * np.exp( -0.5 * (mapto - mapfrom)**2 / sigma_uv[i]**2 ) * grid_binsize[j]
            
            if normalization_check: # redundant normalization
                sum_raw_output = np.sum(raw_output)
                if np.abs(sum_raw_output - 1) > 1e-2:
                    if verbose: print("Warning: the convolution is not perfectly normalized", sum_raw_output, 'at', input_grid[i])
                if sum_raw_output > 0:
                    raw_output = raw_output/sum_raw_output
                else:
                    raw_output = np.zeros(len(input_grid))

            output_weight += raw_output * (input_weight[i] * grid_binsize[i])
        return output_weight/grid_binsize
        
    def QHII_dz(self, y, z, nion, nH, fesc, alpha_B, CHII,H):
        #return 1/nH*scipy.integrate.quad(lambda z: fesc*nion(z)*cm_Mpc**3/((1+z)*H(1/(1+z))*km_Mpc) \
        #*np.exp(-alpha_B*nH*CHII*xe*scipy.integrate.quad(lambda zp: (1+zp)**2/(H(1/(1+zp))*km_Mpc), z0, z)[0]), z0, 50)[0]
        # function that returns dy/dt
        xe = y*(1+YHe/4)
        QHIIdt = fesc*nion/2.938e+73/nH - CHII*alpha_B*nH*(1+z)**3*xe
        QHIIdz = QHIIdt*(-(1+z)*H*km_Mpc)**(-1)
        return QHIIdz


    def QHII(self, rhoM, model, model_H, model_SFR, par1, par2, pool_cpu, f0=None):
        Nion = 10**53.14
        nH = (1-YHe)*Omegab0*(H0/100)**2*1.88e-29/(mP*1000)
        z_int = np.linspace(20, 0, 64)
        #z_int = np.hstack(([20,15],z_int))
        cosmological_library = cosmological_functions(
            1/(1+z_int), model, model_H, par1, par2)
        H = cosmological_library.H_f(1/(1+z_int), model_H, par1, par2)
        H = scipy.interpolate.interp1d(z_int, H, fill_value='extrapolate')
        
        Pk_arr = []
        for i, z_i in enumerate(z_int):
            HMF_library = HMF(1/(1+z_i), model, model_H, par1, par2, 1e8)
            Pk_arr.append(np.array(HMF_library.Pk(1/(1+z_i), model, par1, par2))*h**3)
        k = kvec/h
        print('finished Pk')
        iterable = [(1/(1+z), rhoM, model, model_H,
                          model_SFR, par1, par2, k, Pk_arr[i], f0) for i,z in enumerate(z_int)]
        nion = pool_cpu.starmap(self.n_ion,tqdm(iterable, total=len(z_int)))
        
        if model_SFR == 'double_power':
            sigma_SFRD = 0.1
        elif model_SFR == 'Puebla':
            sigma_SFRD = 0.2
        
        nion = self.convolve_on_grid(z_int,nion,sigma_SFRD)
        nion = Nion*nion
        print('finished nion')
        nion = scipy.interpolate.interp1d(
            z_int, nion, fill_value='extrapolate')
        #QHIIdz = self.QHII_dz(z_int, nion(z_int), nH, fesc, alpha_B, CHII,H(z_int))
        def model(y,t):
            return self.QHII_dz(y,t,nion(t), nH, fesc, alpha_B, CHII,H(t))
        
        QHII = odeint(model,1e-13,z_int)
        return z_int, QHII


    def tau_reio(self, z_span, qhii, model, model_H, par1, par2):
        #z_span, qhii = self.QHII(rhoM, model, model_H, model_SFR, par1, par2, f0=f0)
        a_int = 1/(1+z_span)

        #z_span = np.linspace(0,50,50)
        #a_int = np.linspace(1/51,1,1000)
        cosmological_library = cosmological_functions(
            a_int, model, model_H, par1, par2)
        H = cosmological_library.H_f(a_int, model_H, par1, par2)
        H = scipy.interpolate.interp1d(a_int, H, fill_value='extrapolate')
        #QHII = np.array(self.QHII(1/(1+z_span), rhoM, model, model_H, model_SFR, par1, par2, f0))
        qhii[qhii > 1] = 1
        xe = scipy.interpolate.interp1d(z_span, (1+YHe/4)*qhii, fill_value='extrapolate')
        nH = (1-YHe)*Omegab0*(H0/100)**2*1.88e-29/(mP*1000)
        tau_reio = nH*sigma_T*scipy.integrate.quad(lambda z_span: c*1e5*xe(z_span)*(1+z_span)**2/(H(1/(1+z_span))*km_Mpc),min(z_span), max(z_span))[0]
        return tau_reio
