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
        Hprime = a*dH
        Rprime = a*dRda

        ddRda = (-Hprime/H*Rprime +
                 (1+Hprime/H)*R - Omegam0*a**(-3)*H0**2 /
                 (2*H**2) * mu*(R+a/ai)*delta_nl_a - a*dRda)/a**2
        return [dRda, ddRda]

    def radius_solve(self, model, model_H, par1, par2, deltai, delta_nl):
        Hi = cosmological_library.H_f(ai, model_H, par1, par2)
        # ddeltai = Hi*deltai
        R_arr = scipy.integrate.odeint(self.radius_evolution, [0, -ai*Hi*deltai/(3*(1+deltai))], a_arr, args=(
            model, model_H, par1, par2, deltai, delta_nl), tfirst=False)[:, 0]

        return a_arr, R_arr + a_arr/ai

    def virial_theorem(self, model, model_H, par1, par2, deltai, delta_nl, ac):
        G = 1/(8*np.pi)
        a_arr, R_arr = self.radius_solve(
            model, model_H, par1, par2, deltai, delta_nl)
        H_arr = cosmological_library.H_f(a_arr, model_H, par1, par2)
        dH_arr = cosmological_library.dH_f(a_arr, model_H, par1, par2)
        mu_arr = cosmological_library.mu(a_arr, model, model_H, par1, par2)

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
        U = -3/5*G*mu_arr*M*deltaM/R_arr+3/5*(H_dot+H_arr**2)*M*R_arr**2

        radius = scipy.interpolate.interp1d(
            a_arr, R_arr, fill_value='extrapolate')
        virial = scipy.interpolate.interp1d(
            a_arr, T+1/2*U, fill_value='extrapolate')

        a_turn = a_arr[R_arr.argmax()]
        a_vir = scipy.optimize.fsolve(
            lambda x: virial(x), np.array(a_turn))
        a_vir_predict = scipy.optimize.fsolve(
            lambda x: radius(x)/radius(a_turn)-1/2, ai)
        # , a_arr,T+1/2*U
        return a_turn, a_vir, a_vir_predict, radius(a_vir), a_arr, R_arr
