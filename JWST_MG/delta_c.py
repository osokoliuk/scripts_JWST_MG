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

    """
    Set up the non-linear ODE function.

    Parameters:
        a (float): The scale factor.
        y (list): The state vector containing delta and ddeltada.
        model (str): The cosmological model.
        model_H (str): The Hubble model.
        par1 (float): The first parameter.
        par2 (float): The second parameter.

    Returns:
        list: First and second order derivatives of delta.
    """

    def delta_nl_ODE(self, a, y, model, model_H, par1, par2):
        delta, ddeltada = y
        cosmological_library = cosmological_functions(
            a, model, model_H, par1, par2)
        H = cosmological_library.H_f(a, model_H, par1, par2)
        dH = cosmological_library.dH_f(a, model_H, par1, par2)
        if model == 'LCDM' or model == 'wCDM':
            mu = cosmological_library.mu(a, model, model_H, par1, par2)
            dddeltada = -(3/a+dH/H)*ddeltada + (3*Omegam0*mu)/(2*a**5*H **
                                                               2/H0**2)*delta*(1+delta) + 4*ddeltada**2/(3*(1+delta))
        elif model == 'nDGP':
            Hdot = a*H*dH
            beta = 1 + 2*H*par1/c*(1+Hdot/(3*H**2))
            epsilon = 8/(9*beta**2)*(H0*par1/c)**2*Omegam0*a**(-3)
            RRV = (epsilon*delta)**(-1/3)
            mu = cosmological_library.mu(
                a, model, model_H, par1, par2, type='nonlinear', x=RRV)
            dddeltada = -(3/a+dH/H)*ddeltada + (3*Omegam0*mu)/(2*a**5*H **
                                                               2/H0**2)*delta*(1+delta) + 4*ddeltada**2/(3*(1+delta))
        elif model == 'kmoufl':
            ddeltada = 1
        else:
            raise Exception("Incorrect model specified")
        return [ddeltada, dddeltada]

    """
    Set up the linear ODE function.

    Parameters:
        a (float): The scale factor.
        y (list): The state vector containing delta and ddeltada.
        model (str): The cosmological model.
        model_H (str): The Hubble model.
        par1 (float): The first parameter.
        par2 (float): The second parameter.

    Returns:
        list: First and second order derivatives of delta.
    """
    def delta_l_ODE(self, a, y, model, model_H, par1, par2):
        delta, ddeltada = y
        cosmological_library = cosmological_functions(
            a, model, model_H, par1, par2)
        H = cosmological_library.H_f(a, model_H, par1, par2)
        dH = cosmological_library.dH_f(a, model_H, par1, par2)
        if model == 'LCDM' or model == 'wCDM':
            mu = cosmological_library.mu(a, model, model_H, par1, par2)
            dddeltada = -(3/a+dH/H)*ddeltada + (3*Omegam0*mu) / \
                (2*a**5*H**2/H0**2)*delta
        elif model == 'nDGP':
            mu = cosmological_library.mu(
                a, model, model_H, par1, par2, type='linear')
            dddeltada = -(3/a+dH/H)*ddeltada + (3*Omegam0*mu) / \
                (2*a**5*H ** 2/H0**2)*delta
        elif model == 'kmoufl':
            ddeltada = 1
        else:
            raise Exception("Incorrect model specified")
        return [ddeltada, dddeltada]

    """
    Spherical collapse system given the initial deltai, model, model_H, par1, and par2.
    
    Args:
        deltai (float): The initial delta value.
        model (object): The model object.
        model_H (object): The model_H object.
        par1 (float): The par1 value.
        par2 (float): The par2 value.
    
    Returns:
        float: The value of scale factor at which spherical collapse occurs.
    """

    def collapse(self, deltai, model, model_H, par1, par2):
        dt = 0.00002
        ddeltai = deltai/ai
        init = [deltai, ddeltai]
        system = ode(self.delta_nl_ODE)
        system.set_f_params(*[model, model_H, par1, par2])
        system.set_initial_value(init, ai)

        data = []
        while system.successful() and system.y[0] <= 1e7 and system.t <= 1/(1-0.5):
            data.append([system.t + dt, system.integrate(system.t + dt)[0]])
            # plt.scatter(system.t+dt, system.integrate(system.t + dt)[0], c ='tab:blue')

        data = np.array(data)
        return data

    """
    Perform a binary search to find the index of the element in delta_ini that is closest to the given ac value within a specified absolute error. 

    Args:
        ac (float): The target value for the search.
        model (object): The model object used for prediction.
        model_H (object): The model object used for prediction.
        par1 (float): The parameter for the collapse function.
        par2 (float): The parameter for the collapse function.
        low (int): The lower bound of the search range.
        high (int): The upper bound of the search range.
        abs_err (float): The acceptable absolute error for the search.

    Returns:
        int: The index of the element in delta_ini that is closest to the given ac value within the specified absolute error.
    """

    def binary_search_di(self, ac, model, model_H, par1, par2, low, high, abs_err):
        mid = (low + high)//2
        ac_predict = self.collapse(
            delta_ini[mid], model, model_H, par1, par2)[-1, 0]
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

    """
    Generate the data points for the linear delta based on the given parameters.

    Parameters:
        deltai_collapse (float): The ICS for delta.
        a (float): The upper limit for scale factor.
        model (object): The model object.
        model_H (object): The H model object.
        par1 (float): Parameter 1.
        par2 (float): Parameter 2.

    Returns:
        numpy.ndarray: An array containing the data points.
    """

    def linear(self, deltai_collapse, a, model, model_H, par1, par2):
        dt = 0.00001
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
        if model == 'LCDM' or model == 'wCDM':
            mu = cosmological_library.mu(a, model, model_H, par1, par2)
            dddeltada = -(3/a+dH/H)*ddeltada + (3*Omegam0*mu)/(2*a**5*H **
                                                           2/H0**2)*delta*(1+delta) + 4*ddeltada**2/(3*(1+delta))
        elif model == 'nDGP':
            Hdot = a*H*dH
            beta = 1 + 2*H*par1/c*(1+Hdot/(3*H**2))
            epsilon = 8/(9*beta**2)*(H0*par1/c)**2*Omegam0*a**(-3)
            RRV = (epsilon*delta)**(-1/3)
            mu = cosmological_library.mu(
                a, model, model_H, par1, par2, type='nonlinear', x=RRV)
            dddeltada = -(3/a+dH/H)*ddeltada + (3*Omegam0*mu)/(2*a**5*H **
                                                               2/H0**2)*delta*(1+delta) + 4*ddeltada**2/(3*(1+delta))
        elif model == 'kmoufl':
            ddeltada = 1
        else:
            raise Exception("Incorrect model specified")

        return [ddeltada, dddeltada]

    """
    This function integrates a non-linear ordinary differential equation system using scipy's odeint function.

    Parameters:
        deltai_collapse (float): The collapse parameter
        a (array): Array of values
        model (str): The model used
        model_H (str): The Hubble model
        par1 (int): Parameter 1
        par2 (int): Parameter 2

    Returns:
        array: The first column of the integrated array
    """

    def non_linear(self, deltai_collapse, a, model, model_H, par1, par2):
        ai = a[0]
        ddeltai = deltai_collapse/ai
        init = [deltai_collapse, ddeltai]

        delta_arr = scipy.integrate.odeint(self.delta_nl_ODE2, init, a, args=(
            model, model_H, par1, par2), tfirst=False)

        return delta_arr[:, 0]

    def delta_c_at_ac(self, ac, model, model_H, par1, par2):
        return self.linear(self.binary_search_di(ac, model, model_H, par1, par2, 0, len(delta_ini)-1, abs_err), ac, model, model_H, par1, par2)[-1, 1]
