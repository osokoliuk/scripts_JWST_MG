
from JWST_MG.UVLF import UVLF
from JWST_MG.constants import *

model = 'nDGP'
model_H = 'nDGP'
model_SFR = 'phenomenological_extreme'

par1 = 150000
par2 = 1

Mh0 = 1e12
f0 = 0.05

######################################################################


z = 10
Masses = np.logspace(6.5,17,100)
a = 1/(1+z)
SMD_obs = UVLF(a, model, model_H, model_SFR, par1, par2, Masses, f0)
SFR = SMD_obs.compute_uv_luminosity_function(a, rhom, model, model_H, model_SFR, par1, par2, Masses, f0)
plt.plot(z, SFR) 


my_xlims=np.r_[1e8, 1e12]
my_ylims=np.r_[1e2, 1e8]
#plt.xlim(*my_xlims)
#plt.ylim(*my_ylims)
plt.yscale('log')

plt.savefig('HMF.pdf')
#delta_c_at_ac(self, ac, model, model_H, par1, par2):
