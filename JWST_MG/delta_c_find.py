import cosmological_functions

def delta_nl_ODE(a, y, par1, par2):
    delta,ddeltada = y
    H = H_f(a)
    dH = dH_f(a)
    dddeltada          = -(3/a+dH/H)*ddeltada + (3*Omegam0*mu(par1, par2, a))/(2*a**5*H**2/H0**2)*delta*(1+delta) + 4*ddeltada**2/(3*(1+delta))
    return [ddeltada, dddeltada]

def delta_l_ODE(a, y, par1, par2):
    delta,ddeltada = y
    H = H_f(a)
    dH = dH_f(a)
    dddeltada          = -(3/a+dH/H)*ddeltada + (3*Omegam0*mu(par1, par2, a))/(2*a**5*H**2/H0**2)*delta
    return [ddeltada, dddeltada]

def collapse(deltai, par1, par2):
    ai = 1e-5
    dt = 0.0001
    ddeltai = deltai/ai
    init = [deltai, ddeltai]
    system = ode(delta_nl_ODE)
    system.set_f_params(*[par1, par2])
    system.set_initial_value(init, ai)
    
    data = [0,0]
    while system.successful() and system.y[0] <= 1e7 and system.t <= 1/(1-0.5):
        data = [system.t + dt, system.integrate(system.t + dt)[0]]   
        #plt.scatter(system.t+dt, system.integrate(system.t + dt)[0], c ='tab:blue')
    data = np.array(data)
    ac = data[0]

    return ac

def interpolate_ac(ac, par1, par2):
    delta_i = np.logspace(-4,-1,1000)
    arr_di_ac = []
    for i in tqdm(range(len(delta_i))):
        arr_di_ac.append([delta_i[i],collapse(delta_i[i],pars1,pars2)])
    interp_di_ac = scipy.interpolate.interp1d(np.array(arr_di_ac)[:,1],np.array(arr_di_ac)[:,0],fill_value = 'extrapolate')
    return interp_di_ac(ac)

def linear(deltai_collapse, a, par1, par2):
    ai = 1e-5
    dt = 0.0001
    ddeltai = deltai_collapse/ai    
    init = [deltai_collapse, ddeltai]
    system = ode(delta_l_ODE)
    system.set_f_params(*[par1, par2])
    system.set_initial_value(init, ai)
    
    data = []
    while system.successful() and system.t <= a:
        data.append([system.t + dt, system.integrate(system.t + dt)[0]])   
        #plt.scatter(system.t+dt, system.integrate(system.t + dt)[0], c ='tab:orange')
    
    data = np.array(data)
    return data

def delta_c_at_ac(ac, par1, par2):
    return linear(interp_di_ac(ac),ac,par1,par2)[-1,1]
