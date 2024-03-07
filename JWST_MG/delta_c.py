from JWST_MG.constants import *
from JWST_MG.cosmological_functions import cosmological_functions


class delta_c:
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

    def delta_nl_ODE(self, a, y, model, model_H, par1, par2):
        delta, ddeltada = y
        cosmological_library = cosmological_functions(
            a, model, model_H, par1, par2)
        H = cosmological_library.H_f(a, model_H, par1, par2)
        dH = cosmological_library.dH_f(a, model_H, par1, par2)
        mu = cosmological_library.mu(a, model, model_H, par1, par2)
        dddeltada = -(3/a+dH/H)*ddeltada + (3*Omegam0*mu)/(2*a**5*H **
                                                           2/H0**2)*delta*(1+delta) + 4*ddeltada**2/(3*(1+delta))
        return [ddeltada, dddeltada]

    def delta_l_ODE(self, a, y, model, model_H, par1, par2):
        delta, ddeltada = y
        cosmological_library = cosmological_functions(
            a, model, model_H, par1, par2)
        H = cosmological_library.H_f(a, model_H, par1, par2)
        dH = cosmological_library.dH_f(a, model_H, par1, par2)
        mu = cosmological_library.mu(a, model, model_H, par1, par2)
        dddeltada = -(3/a+dH/H)*ddeltada + (3*Omegam0*mu) / \
            (2*a**5*H**2/H0**2)*delta
        return [ddeltada, dddeltada]

    def collapse(self, deltai, model, model_H, par1, par2):
        dt = 0.00005
        ddeltai = deltai/ai
        init = [deltai, ddeltai]
        system = ode(self.delta_nl_ODE)
        system.set_f_params(*[model, model_H, par1, par2])
        system.set_initial_value(init, ai)

        data = [0, 0]
        while system.successful() and system.y[0] <= 1e7 and system.t <= 1/(1-0.5):
            data = [system.t + dt, system.integrate(system.t + dt)[0]]
            # plt.scatter(system.t+dt, system.integrate(system.t + dt)[0], c ='tab:blue')
        data = np.array(data)
        ac = data[0]

        return ac

    def binary_search_di(self, ac, model, model_H, par1, par2, low, high, abs_err):
        mid = (low + high)//2
        ac_predict = self.collapse(delta_ini[mid], model, model_H, par1, par2)
        print(ac_predict)
        if high >= low:
            if abs(ac_predict-ac)/ac_predict <= abs_err:
                return delta_ini[mid]
            elif ac_predict < ac:
                return self.binary_search_di(ac, model, model_H, par1, par2, low, mid - 1, abs_err)
            else:
                return self.binary_search_di(ac, model, model_H, par1, par2, mid + 1, high, abs_err)
        else:
            low = 0
            high = len(delta_ini) - 1
            return self.binary_search_di(ac, model, model_H, par1, par2, low, high, abs_err*2)


    def linear(self, deltai_collapse, a, model, model_H, par1, par2):
        dt = 0.000001
        ddeltai = deltai_collapse/ai
        init = [deltai_collapse, ddeltai]
        system = ode(self.delta_l_ODE)
        system.set_f_params(*[model, model_H, par1, par2])
        system.set_initial_value(init, ai)

        data = []
        while system.successful() and system.t <= a:
            data.append([system.t + dt, system.integrate(system.t + dt)[0]])
            # plt.scatter(system.t+dt, system.integrate(system.t + dt)[0], c ='tab:orange')

        data = np.array(data)
        return data

    def delta_nl_ODE2(self, y, a, model, model_H, par1, par2):
        delta, ddeltada = y
        cosmological_library = cosmological_functions(
            a, model, model_H, par1, par2)
        H = cosmological_library.H_f(a, model_H, par1, par2)
        dH = cosmological_library.dH_f(a, model_H, par1, par2)
        mu = cosmological_library.mu(a, model, model_H, par1, par2)
        dddeltada = -(3/a+dH/H)*ddeltada + (3*Omegam0*mu)/(2*a**5*H **
                                                           2/H0**2)*delta*(1+delta) + 4*ddeltada**2/(3*(1+delta))
        return [ddeltada, dddeltada]

    def non_linear(self, deltai_collapse, a, model, model_H, par1, par2):
        ai = a[0]
        ddeltai = deltai_collapse/ai
        init = [deltai_collapse, ddeltai]

        delta_arr = scipy.integrate.odeint(self.delta_nl_ODE2, init, a, args=(
            model, model_H, par1, par2), tfirst=False)

        return delta_arr[:, 0]


    def delta_c_at_ac(self, ac, model, model_H, par1, par2):
        return self.linear(self.binary_search_di(ac, model, model_H, par1, par2, 0, len(delta_ini), abs_err), ac, model, model_H, par1, par2)[-1, 1]
