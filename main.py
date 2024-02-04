
from JWST_MG.UVLF import UVLF
from JWST_MG.constants import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import commah

model = 'nDGP'
model_H = 'nDGP'
model_SFR = 'phenomenological_extreme'

par1 = 150000
par2 = 1

Mh0 = 1e12
f0 = 0.05

######################################################################


z = np.array([0,1,2,3,4,5,6,7,8])
Masses = 1e8
a = 1/(1+z)
SMD_obs = UVLF(a, model, model_H, model_SFR, par1, par2, Masses, f0)
Mh_EPS, dMh_EPS = SMD_obs.Mh_EPS(a, rhom, model, model_H, par1, par2, Mh0)
plt.plot(np.log10(1+z), np.log10(dMh_EPS)) 



def halo_accretion_rate(mhalo, redshift):
    # Fakhouri 2010
    # Mhalo in Msun
    mhalo_dot = 46.1 * (1 + 1.11*redshift) * np.sqrt(Omegam0*(1+redshift)**3 + (1-Omegam0))  \
     * (mhalo / 1e12)**(1.1)
    #mhalo_dot = 25.3 * (1 + 1.65*redshift) * np.sqrt(Omegam0*(1+redshift)**3 + (1-Omegam0))  \
    # * (mhalo / 1e12)**(1.1)
    #corr = 10**(-0.1) # down 0.1 dex consider the drop of sigma8
    return mhalo_dot 



#dMhdt_mean = 25.3*(Mh0/1e12)**1.1*(1+1.65*z)*np.sqrt(Omegam0*(1+z)**3+1-Omegam0)
plt.plot(np.log10(1+z), np.log10(halo_accretion_rate(Mh0, z))) 



my_xlims=np.r_[1e8, 1e12]
my_ylims=np.r_[1e2, 1e8]
#plt.xlim(*my_xlims)
#plt.ylim(*my_ylims)
#plt.yscale('log')

plt.savefig('HMF.pdf')
#delta_c_at_ac(self, ac, model, model_H, par1, par2):
