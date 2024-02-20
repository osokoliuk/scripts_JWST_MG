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
        Hdot = a*H*dH
        rho_m = 3*H0**2*Omegam0*a**(-3)
        delta_nl = deltac_library.non_linear(deltac_library.interpolate_ac(
            ac, model, model_H, par1, par2), a, model, model_H, par1, par2)[-1, -1]
        ddRda = ((Hdot+H**2-mu*rho_m*delta_nl/6) -
                 a*(H**2+Hdot)*dRda/R)/(a**2*H**2)
        return [dRda, ddRda]
