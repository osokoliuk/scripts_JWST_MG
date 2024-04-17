import sys
sys.path.insert(0, "../")
from JWST_MG.constants import *
from JWST_MG.cosmological_functions import cosmological_functions
from JWST_MG.delta_c import delta_c
from JWST_MG.HMF import HMF
from JWST_MG.SMF import SMF
from JWST_MG.SMD import SMD
from JWST_MG.UVLF import UVLF
plt.rcParams.update({"text.usetex": True})



x = [[],[],[],[],[],[],[]]
y = [[],[],[],[],[],[],[]]
yerr = [[],[],[],[],[],[],[]]

zs = [4,5,6,7,8,9,10]

for i in range(len(zs)):
    obs = number_density(feature='GSMF', z_target=zs[i], h=h)
    j_data = 0
    k_func = 0

    for ii in range(obs.n_target_observation):
        data       = obs.target_observation['Data'][ii]
        datatype   = obs.target_observation['DataType'][ii]
        data[:,1:] = data[:,1:]
        if datatype == 'data':
            x[i] = np.concatenate((x[i], 10**data[:,0]), axis=None)
            y[i] = np.concatenate((y[i], data[:,1]), axis=None)
            yerr[i] = np.concatenate((yerr[i], data[:,1]-data[:,3]+data[:,2]- data[:,1]), axis=None)

            j_data +=1




path = '../observational_data/GSMF'
zs2 = [4,5,6,7,8]
for i in range(len(zs2)):
    Navarro = np.loadtxt(path + "/Navarro_z"+str(zs2[i])+".dat")
    x[i] = np.concatenate((x[i], 10**Navarro[:,0]), axis=None)
    y[i] = np.concatenate((y[i], 1e-4*Navarro[:,1]), axis=None)
    yerr[i] = np.concatenate((yerr[i], 2*1e-4*Navarro[:,2]), axis=None)


model = 'nDGP'
model_H = 'nDGP'
model_SFR = 'Behroozi'
par2 = 1
f0 = 0.1
def SMF_func(z, par1):
    SMF_library = SMF(1/(1+z), model, model_H, model_SFR, par1, par2, 1e8)
    HMF_library = HMF(1/(1+z), model, model_H, par1, par2, 1e8)
    Pk = np.array(HMF_library.Pk(1/(1+z), model, par1, par2))*h**3
    k = kvec/h
    Masses = np.logspace(8,14,100)
    Masses_star, SMF_sample = SMF_library.SMF_obs(Masses, rhom, 1/(1+z), model_H, model, model_SFR, par1, par2, k, Pk, f0)
    return Masses_star, SMF_sample


def SMF_pool(log_par1):
    par1 = 10**log_par1
    for k, zi in enumerate(zs):
        Masses_star, SMF_sample = SMF_func(zi, par1)
        y_th = scipy.interpolate.interp1d(Masses_star, SMF_sample, fill_value='extrapolate')(x[k])
        sigma2 = yerr[k]**2
        result[i] += -0.5 * np.sum((y[k] - y_th) ** 2 / sigma2 + np.log(sigma2))
    return result