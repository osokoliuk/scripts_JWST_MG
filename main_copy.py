
from colossus.cosmology import cosmology
from JWST_MG.UVLF import UVLF
from JWST_MG.HMF import HMF

from JWST_MG.constants import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import commah
ax = plt.subplot(111)
model = 'nDGP'
model_H = 'nDGP'
model_SFR = 'toy'

par1 = 150000
par2 = 1

Mh0 = 1e12
f0 = 0.05


z = 10
Masses = np.logspace(8, 16, 50)
a = 1/(1+z)

x = np.loadtxt("/home/oleksii/highz-empirical-variability/data.txt")[:, 0]
y = np.loadtxt("/home/oleksii/highz-empirical-variability/data.txt")[:, 1]
y = scipy.interpolate.interp1d(np.log10(x), y, fill_value='extrapolate')


HMF_library = HMF(a, model, model_H, par1, par2, Masses)
phi_halo_arr = HMF_library.ST_mass_function(
    rhom, Masses, a, model_H, model, par1, par2)  # mapfunc_mhalo_to_muv(self, a, rhoM, model, model_H, model_SFR, par1, par2, Masses, f0, dust_norm="fixed", include_dust=True):

plt.plot(np.log10(Masses), phi_halo_arr)

cosmology.setCosmology('planck18');
from colossus.lss import mass_function
mfunc = mass_function.massFunction(
    Masses, 10, mdef='fof', model='sheth99', q_out='dndlnM', sigma_args={'filt': 'tophat'})
Masses = Masses*h
mfunc = mfunc/h**1/Masses

plt.plot(np.log10(Masses), mfunc)

plt.plot(np.log10(Masses), y(np.log10(Masses)), ls=':', c='tab:green', lw=3)
plt.yscale('log')
plt.savefig('HMF.pdf', bbox_inches='tight')
