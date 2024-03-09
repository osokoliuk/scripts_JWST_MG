from JWST_MG.constants import *
from JWST_MG.cosmological_functions import cosmological_functions
from JWST_MG.delta_c import delta_c

class reionization:
    ########################################################################
    # Initialize a class delta_c (critical threshold for linear density perturbations)
    # float ac - scale factor, at which spherical collapse happens
    # string model - model of MG for the derivation of mu parameter
    # string model_H - model of MG for H(a)
    # float par1, par2 - corresponding MG parameters
    # float delta_i - initial linear/non-linear overdensity at a = 1e-5
    ########################################################################

    def __init__(self, a_arr, model, model_H, par1, par2):
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

        delta_nl = self.delta_nl_a(a)
        H = cosmological_library.H_f(a, model_H, par1, par2)
        dH = cosmological_library.dH_f(a, model_H, par1, par2)
        mu = cosmological_library.mu(a, model, model_H, par1, par2)
        Hprime = a*dH
        Rprime = a*dRda

        if model_H == "LCDM" or model_H == "wCDM":
            ddRda = (-Hprime/H*Rprime +
                 (1+Hprime/H)*R - Omegam0*a**(-3)*H0**2 /
                 (2*H**2) * mu*(R+a/ai)*delta_nl - a*dRda)/a**2
        elif model_H == "nDGP" or "kmoufl":
            ddRda = 1
        else:
            raise ValueError("Unknown cosmology type!")

        return [dRda, ddRda]

    def radius_solve(self, model, model_H, par1, par2, a_arr):
        deltai = self.deltai
        Hi = cosmological_library.H_f(ai, model_H, par1, par2)
        R_arr = scipy.integrate.odeint(self.radius_evolution, [0, -ai*Hi*deltai/(3*(1+deltai))], a_arr, args=(
            model, model_H, par1, par2, a_arr), tfirst=False)[:, 0]

        return R_arr + a_arr/ai

    def virial_theorem(self, model, model_H, par1, par2, a_arr):
        G = 1/(8*np.pi)
        ac = a_arr[-1]
        R_arr = self.radius_solve(
            model, model_H, par1, par2, a_arr)
        H_arr = cosmological_library.H_f(a_arr, model_H, par1, par2)
        dH_arr = cosmological_library.dH_f(a_arr, model_H, par1, par2)
        mu_arr = cosmological_library.mu(a_arr, model, model_H, par1, par2)
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

        if model_H == "LCDM":
            U = -3/5*G*mu_arr*M*deltaM/R_arr+3/5*(H_dot+H_arr**2)*M*R_arr**2
            virial = T+1/2*U
        elif model_H == "wCDM":
            rho_Lambda = 3*H0**2*(1-Omegam0-Omegar0)*a_arr**(3*(1+par1))
            U_Lambda = -4/5*np.pi*G*rho_Lambda*M*R**2
            U_N = -3/5*G*M**2/R
            virial = T+1/2*U_N - U_Lambda
        elif model_H == "nDGP" or model_H == "kmoufl":
            U = 1
            virial = U
        else:
            raise ValueError("Unknown cosmology type!")

        virial[virial < 0] = 1e60

        virial = scipy.interpolate.interp1d(
            a_arr, virial, fill_value='extrapolate')

        a_turn = a_arr[R_arr.argmax()]
        a_vir = scipy.optimize.minimize_scalar(
            lambda x: virial(x), bounds=(ai, self.ac), method="bounded").x

        return a_turn, a_vir

    def Delta_vir(self, model, model_H, par1, par2, a_arr):
        ac = a_arr[-1]

        if model_H == "LCDM":
            a_turn, a_vir = self.virial_theorem(
            model, model_H, par1, par2, a_arr)
            Deltavir = (1+self.delta_nl_a(a_vir))*(ac/a_vir)**3
        else:
            raise Exception("Incorrect model specified")
        return a_vir, Deltavir

    def minimum_Mhalo(model, model_H, par1, par2, a_arr):
        Delta_vir = self.Delta_vir(
            model, model_H, par1, par2, a_arr)
        rho_vir = Delta_vir*rhom
        R_vir = (3*GN*M/(4*np.pi*rho_vir))**(1/3)
        def v_vir(Mhalo): return np.sqrt(GN*Mhalo/R_vir)
        Mhalo_min = scipy.optimize.minimize_scalar(
            lambda x: 0.75*mu_mol*mH*v_vir(x)**2/(2*kB) - 4000, bounds=(1e6, 1e16), method="bounded")
        Mhalo_min = Mhalo_min.x

        return Mhalo_min

    # def n_ion():
    #    Nion =
