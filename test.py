from JWST_MG.UVLF import UVLF
from JWST_MG.HMF import HMF
from JWST_MG.reionization import reionization
from JWST_MG.delta_c import delta_c

from JWST_MG.constants import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

model = 'wCDM'
model_H = 'wCDM'
model_SFR = 'toy'

par2 = 6/11
par1 = -0.7

ac_arr = np.linspace(1, 1/(1+2.5), 15)

for ac in ac_arr:
    deltac = delta_c(ac, model, model_H, par1, par2)
    plt.scatter(1/ac-1, deltac.delta_c_at_ac(ac, model,
                model_H, par1, par2), c='tab:blue')


model = 'LCDM'
model_H = 'LCDM'
model_SFR = 'toy'

par2 = 6/11
par1 = -1.3

for ac in ac_arr:
    deltac = delta_c(ac, model, model_H, par1, par2)
    plt.scatter(1/ac-1, deltac.delta_c_at_ac(ac, model,
                model_H, par1, par2), c='tab:orange')


# plt.xscale('log')
plt.savefig('HMF.pdf')
