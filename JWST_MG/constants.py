from scipy.integrate import ode, odeint
import scipy.constants as SPC
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
import numpy as np
from numpy import inf
from classy import Class
from scipy.optimize import fsolve
import math
import scipy
from tqdm import tqdm
import integration_library as IL
from scipy.interpolate import interp1d
import classy
import sys
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import pandas as pd
import scipy.integrate as intg
import astropy.constants as con
import pickle
import os
from astrodatapy.number_density import number_density
from astrodatapy.correlation    import correlation
from astrodatapy.clustering     import clustering
from multiprocessing import Pool
import emcee
from parallelbar import progress_map, progress_starmap
from scipy.stats import qmc
from numpy import log, log10, sqrt, exp, pi, arange, logspace, linspace, r_, c_

os.environ["OMP_NUM_THREADS"] = "1"


kmfl_settings = {'A_s': 2.101e-9,
                 'n_s': 0.9665,
                 'tau_reio': 0.0561,
                 'omega_b': 0.02242,
                 'omega_cdm': 0.11933,
                 'h': 0.6766,
                 'YHe': 0.2425,
                 'T_cmb': 2.7255,
                 'gauge': 'newtonian',  # FOR MGCLASS TO WORK, GAUGE NEEDS TO BE NEWTONIAN
                 'k_pivot': 0.05,
                 'mg_z_init': 99.000,
                 'l_logstep': 1.025,
                 'l_linstep': 15,
                 'P_k_max_1/Mpc': 1000.0,
                 'l_switch_limber': 9,
                 'perturb_sampling_stepsize': 0.05,
                 'output': 'tCl,pCl,lCl,mPk',
                 'l_max_scalars': 3000,
                 'lensing': 'yes',
                 'mg_ansatz': 'kmoufl'}

H0 = 67.66
Omegam0 = (0.02242/(H0/100)**2+0.11933/(H0/100)**2)
Omegar0 = 8.493e-5
Omegab0 = 0.02242/(H0/100)**2


c = 299792.45800000057
Tcmb = 2.72e6
YHe = 0.243
kB = 1.380649*1e-23
mP = 1.6726e-27
Msol = 1.998e30
mP_Msol = mP/Msol
mu_mol = 1.22
h = H0/100
GN = 4.301*10**(-9)
rho = 3*H0**2*Omegam0/(8*np.pi*GN)
rhocr = 2.77536627e11
rhom = rhocr*Omegam0
CHII = 3
THII = 2*1e4
fesc = 0.1
cm_Mpc = 3.24078e-25
km_Mpc = 3.24078e-20
MsolMpc3_to_gcm3 = 6.77e-23
alpha_B = 2.5*1e-13
sigma_T = 6.6525e-25

kvec = np.logspace(np.log10(0.000001), np.log10(1000.0), 10000)
delta_ini = np.logspace(-5, 0, 100000)
abs_err = 1e-5
ai = 1e-5
beta_arr = np.linspace(0, 0.5, 15)
K0_arr = np.linspace(0.1, 1, 15)

H_arr_kmoufl = np.zeros(shape=(15, 15), dtype=object)
dH_arr_kmoufl = np.zeros(shape=(15, 15), dtype=object)
H_int_kmoufl = np.zeros(shape=(15, 15), dtype=object)
dH_int_kmoufl = np.zeros(shape=(15, 15), dtype=object)

M = {}
path = os.path.dirname(os.path.realpath(__file__))
kmoufl_H = np.load(path  + "/kmoufl_H.npy", allow_pickle=True)
kmoufl_dH = np.load(path + "/kmoufl_dH.npy", allow_pickle=True)
