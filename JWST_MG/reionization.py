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

    def __init__(self, ac, model, model_H, par1, par2):
        self.ac = ac
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

    def radius_evolution(self, y, a, model, model_H, par1, par2, deltai, delta_nl):
        R, dRda = y
        delta_nl_a = scipy.interpolate.interp1d(
            a_arr, delta_nl, fill_value='extrapolate')
        delta_nl_a = delta_nl_a(a)
        H = cosmological_library.H_f(a, model_H, par1, par2)
        dH = cosmological_library.dH_f(a, model_H, par1, par2)
        mu = cosmological_library.mu(a, model, model_H, par1, par2)
        ddRda = -dH/H*dRda + (1+dH/H)*R-Omegam0*a**(-3)*H0**2 / (2*H**2) \
            * mu*delta_nl_a*(R+a/ai)
        return [dRda, ddRda]

    def radius_solve(self, model, model_H, par1, par2, deltai, delta_nl):
        # Hi = cosmological_library.H_f(ai, model_H, par1, par2)
        # ddeltai = Hi*deltai
        R_arr = scipy.integrate.odeint(self.radius_evolution, [0, -deltai/3], a_arr, args=(
            model, model_H, par1, par2, deltai, delta_nl), tfirst=False)[:, 0]

        return a_arr, R_arr + a_arr/ai
