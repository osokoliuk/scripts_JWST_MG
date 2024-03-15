from JWST_MG.UVLF import UVLF
from JWST_MG.HMF import HMF
from JWST_MG.reionization import reionization
from JWST_MG.delta_c import delta_c

from JWST_MG.constants import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import commah
ax = plt.subplot(111)
model = 'nDGP'
model_H = 'nDGP'
model_SFR = 'toy'

closest_beta = (np.absolute(K0_arr-0.5)).argmin()

print(closest_beta)
