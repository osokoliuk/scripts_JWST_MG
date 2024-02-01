
from JWST_MG.HMF import HMF
from JWST_MG.constants import *
model = 'nDGP'
model_H = 'nDGP'

par1 = 150000
par2 = 1

HMF_fid = HMF(1, kvec/h, model, model_H, par1, par2)
z = 15
Masses = np.logspace(1e12,1e15,15)
print(HMF_fid.ST_mass_function(kvec/h,rhom, Masses, 1/(1+z), model_H, model, par1, par2))

#delta_c_at_ac(self, ac, model, model_H, par1, par2):
