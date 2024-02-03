
from JWST_MG.SMD import SMD
from JWST_MG.constants import *

model = 'nDGP'
model_H = 'nDGP'
model_SFR = 'phenomenological_extreme'

par1 = 150000
par2 = 1

Masses = np.logspace(7,16,50)
f0 = 0.05

######################################################################

x = np.loadtxt("/home/oleksii/smd_lcdm.txt")[:,0]
y = np.loadtxt("/home/oleksii/smd_lcdm.txt")[:,1]
plt.loglog(x,y)
z_arr = np.array([9.1])
for i in range(len(z_arr)):
    z = z_arr[i]
    a = 1/(1+z)
    SMD_obs = SMD(a, model, model_H, model_SFR, par1, par2, Masses, f0)
    Masses_star, SMD_obs = SMD_obs.SMD(Masses, rhom, a, model_H, model,model_SFR, par1, par2, f0)
    plt.loglog(Masses_star, SMD_obs) 



my_xlims=np.r_[1e8, 1e12]
my_ylims=np.r_[1e2, 1e8]
plt.xlim(*my_xlims)
plt.ylim(*my_ylims)

plt.savefig('HMF.pdf')
#delta_c_at_ac(self, ac, model, model_H, par1, par2):
