
from JWST_MG.HMF import HMF
from JWST_MG.constants import *
model = 'nDGP'
model_H = 'nDGP'

par1 = 150000
par2 = 1

Masses = np.logspace(6.5,15,50)

HMF_fid = HMF(1, kvec/h, model, model_H, par1, par2, Masses)
z = 0
plt.loglog(Masses,HMF_fid.ST_mass_function(kvec/h,rhom, Masses, 1/(1+z), model_H, model, par1, par2))
plt.savefig('HMF.pdf')
#delta_c_at_ac(self, ac, model, model_H, par1, par2):
