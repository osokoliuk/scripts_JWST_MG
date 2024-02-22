from JWST_MG.constants import *
from JWST_MG.cosmological_functions import cosmological_functions

class reionization:
    ########################################################################
    # Initialize a class delta_c (critical threshold for linear density perturbations)
    # float ac - scale factor, at which spherical collapse happens
    # string model - model of MG for the derivation of mu parameter
    # string model_H - model of MG for H(a)
    # float par1, par2 - corresponding MG parameters
    # float delta_i - initial linear/non-linear overdensity at a = 1e-5
    ########################################################################

    def __init__(self, ac, model, model_H, par1, par2):
        self.ac = ac
        self.model = model
        self.model_H = model_H
        self.par1 = par1
        self.par2 = par2

    def radius_evolution(self, a, model, model_H, par1, par2):
        R, dRda = y
        cosmological_library = cosmological_functions(
            a, model, model_H, par1, par2)
        deltac_library = cosmological_functions(
            a, model, model_H, par1, par2)
        H = cosmological_library.H_f(a, model_H, par1, par2)
        dH = cosmological_library.dH_f(a, model_H, par1, par2)
        mu = cosmological_library.mu(a, model, model_H, par1, par2)
        deltai = deltac_library.interpolate_ac(
            ac, model, model_H, par1, par2)
        delta_nl = deltac_library.non_linear(
            deltai, a, model, model_H, par1, par2)[-1, -1]
        Hdot = a*H*dH
        G = 1/8*np.pi
        rho_m = 3*H0**2*Omegam0*a**(-3)
        ddRda = ((Hdot+H**2-4*np.pi*G*mu*rho_m, *delta_nl/3) -
                 a*(H**2+Hdot)*dRda)/(a**2*H**2*R)
        return [dRda, ddRda]

    def radius_solve(self, model, model_H, par1, par2):
        a_arr = np.linspace(1, 1e-5, 1000)
        R_arr = scipy.integrate.odeint(radius_evolution, [1e-7, 1], a_arr, args=(
            model, model_H, par1, par2,), tfirst=False)[:, 0]

        return a_arr, R_arr
