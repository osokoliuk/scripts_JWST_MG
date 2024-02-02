
from JWST_MG.SMF import SMF
from JWST_MG.constants import *
model = 'nDGP'
model_H = 'nDGP'
model_SFR = 'double_power'

par1 = 150000
par2 = 1

Masses = np.logspace(6.5,15,50)
f0 = 0.05

SMF_fid = SMF(1, model, model_H, model_SFR, par1, par2, Masses, f0)

z = 0
Masses_star, SMF_obs = SMF_fid.SMF_obs(Masses, rhom, 1, model_H, model,model_SFR, par1, par2, f0)
plt.loglog(Masses_star, SMF_obs)
plt.savefig('HMF.pdf')
#delta_c_at_ac(self, ac, model, model_H, par1, par2):
a