
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


z = np.array([0,1,2,3,4,5,6,7,8])
Masses = 1e8
a = 1/(1+z)
SMD_obs = UVLF(a, model, model_H, model_SFR, par1, par2, Masses, f0)
Mh, dMhdt = SMD_obs.Mh_EPS(a, rhom, model_H, model, par1, par2, Mh0)
plt.plot(z, dMhdt) 



my_xlims=np.r_[1e8, 1e12]
my_ylims=np.r_[1e2, 1e8]
#plt.xlim(*my_xlims)
#plt.ylim(*my_ylims)
plt.yscale('log')

plt.savefig('HMF.pdf')
#delta_c_at_ac(self, ac, model, model_H, par1, par2):
