from scipy.integrate import ode, odeint
import scipy.constants as SPC
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
import math
import scipy
from tqdm import tqdm 

common_settings = {'A_s':2.101e-9,
        'n_s':0.9665,
        'tau_reio':0.0561,
        'omega_b':0.02242,
        'omega_cdm':0.11933,
        'h':0.6766,
        'YHe':0.2425,
        'T_cmb':2.7255,
        'gauge':'newtonian', #FOR MGCLASS TO WORK, GAUGE NEEDS TO BE NEWTONIAN
        'k_pivot': 0.05,
        'mg_z_init': 111111.000,
        'l_logstep': 1.025,
        'l_linstep':15,
        'P_k_max_1/Mpc':3.0,
        'l_switch_limber':9,
        'perturb_sampling_stepsize': 0.05,
        'output':'tCl,pCl,lCl,mPk',
        'l_max_scalars': 3000,
        'lensing': 'yes',
        'mg_ansatz':'kmoufl'}


H0 = 67.66
Omegam0 = (0.02242/(H0/100)**2+0.11933/(H0/100)**2)
Omegar0 = 8.493e-5
c = 299792.45800000057


K0 = [0, 0.5, 1]
beta = np.linspace(0,0.5,10)
H_arr_kmoufl = [[],[],[]]
dH_arr_kmoufl = [[],[],[]]
H_int_kmoufl = [[],[],[]]
dH_int_kmoufl = [[],[],[]]

M = {}