import sys
sys.path.insert(0, "../")
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import odeint
from scipy import integrate
from scipy.special import lambertw
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.pylab as pl
import matplotlib as mpl
import matplotlib
from scipy import interpolate
import scipy
from JWST_MG.UVLF import UVLF
from JWST_MG.HMF import HMF
from JWST_MG.reionization import reionization
from JWST_MG.delta_c import delta_c

from JWST_MG.constants import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

for i in range(len(K0_arr)):
    for j in tqdm(range(len(beta_arr))):
        kmfl_settings['beta_kmfl'] = beta_arr[j]
        kmfl_settings['k0_kmfl'] = K0_arr[i]

        cosmo_kmfl = Class()
        cosmo_kmfl.set(kmfl_settings)
        cosmo_kmfl.compute()
        a = np.logspace(-6, 0, 10000)
        H_arr_kmoufl[i, j] = ([cosmo_kmfl.Hubble(1/ai-1)*c for ai in a])
        dH_arr_kmoufl[i, j] = (np.gradient(H_arr_kmoufl[i][j])/np.gradient(a))
        H_int_kmoufl[i, j] = (scipy.interpolate.interp1d(
            a, H_arr_kmoufl[i][j], fill_value='extrapolate'))
        dH_int_kmoufl[i, j] = (scipy.interpolate.interp1d(
            a, dH_arr_kmoufl[i][j], fill_value='extrapolate'))

np.save('kmoufl_H', H_int_kmoufl)
np.save('kmoufl_dH', dH_int_kmoufl)
