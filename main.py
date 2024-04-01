from JWST_MG.UVLF import UVLF
from JWST_MG.HMF import HMF
from JWST_MG.reionization import reionization
from JWST_MG.delta_c import delta_c

from JWST_MG.constants import *
from JWST_MG.cosmological_functions import cosmological_functions

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import commah
ax = plt.subplot(111)
model = 'kmoufl'
model_H = 'kmoufl'
model_SFR = 'toy'


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
                 'mg_z_init': 111111.000,
                 'l_logstep': 1.025,
                 'l_linstep': 15,
                 'P_k_max_1/Mpc': 3.0,
                 'l_switch_limber': 9,
                 'perturb_sampling_stepsize': 0.05,
                 'output': 'tCl,pCl,lCl,mPk',
                 'l_max_scalars': 3000,
                 'lensing': 'yes',
                 'mg_ansatz': 'nDGP'}

kmfl_settings['rc'] = 1e8
# kmfl_settings['beta_kmfl'] = 0.1
# kmfl_settings['k0_kmfl'] = 0.5


plt.savefig('delta_c.pdf')
