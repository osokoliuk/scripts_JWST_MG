import sys
sys.path.insert(0, "../")
from JWST_MG.constants import *
from JWST_MG.cosmological_functions import cosmological_functions
from JWST_MG.delta_c import delta_c
from JWST_MG.HMF import HMF
from JWST_MG.SMF import SMF
from JWST_MG.SMD import SMD
from JWST_MG.UVLF import UVLF
import zeus
plt.rcParams.update({"text.usetex": True})
from tqdm.contrib.concurrent import process_map  # or thread_map


x = [[],[],[],[],[],[],[],[],[],[],[]]
y = [[],[],[],[],[],[],[],[],[],[],[]]
yerr = [[],[],[],[],[],[],[],[],[],[],[]]

zs = [0,1,1.75,4,5,6,7,8]

for i in range(len(zs)):
    obs = number_density(feature='GSMF', z_target=zs[i], h=h)
    j_data = 0
    k_func = 0
    if zs[i] == 8:
        pass
    else:
        for ii in range(obs.n_target_observation):
            data       = obs.target_observation['Data'][ii]
            datatype   = obs.target_observation['DataType'][ii]
            if datatype == 'data':
                x[i] = np.concatenate((x[i], 10**data[:,0]), axis=None)
                y[i] = np.concatenate((y[i], data[:,1]), axis=None)
                yerr[i] = np.concatenate((yerr[i], data[:,1]-data[:,3]+data[:,2]- data[:,1]), axis=None)

                j_data +=1
    
    if zs[i] in [4,5,6,7,8]:
        path = '../observational_data/GSMF'
        Navarro = np.loadtxt(path + "/Navarro_z"+str(zs[i])+".dat")
        x[i] = np.concatenate((x[i], 10**Navarro[:,0]), axis=None)
        y[i] = np.concatenate((y[i], 1e-4*Navarro[:,1]), axis=None)
        yerr[i] = np.concatenate((yerr[i], 2*1e-4*Navarro[:,2]), axis=None)


model = 'nDGP'
model_H = 'nDGP'
model_SFR = 'Puebla'
par2 = 1
f0 = 0.1



def SMF_func(z, par1):
    SMF_library = SMF(1/(1+z), model, model_H, model_SFR, par1, par2, 1e8)
    HMF_library = HMF(1/(1+z), model, model_H, par1, par2, 1e8)
    Pk = np.array(HMF_library.Pk(1/(1+z), model, par1, par2))*h**3
    k = kvec/h
    Masses = np.logspace(6,16,100)
    Masses_star, SMF_sample = SMF_library.SMF_obs(Masses, rhom, 1/(1+z), model_H, model, model_SFR, par1, par2, k, Pk, f0)
    return Masses_star, SMF_sample


def SMF_single(log_par1):
    result = 0
    par1 = 10**log_par1
    for k, zi in enumerate(zs):
        Masses_star, SMF_sample = SMF_func(zi, par1)
        y_th = scipy.interpolate.interp1d(Masses_star, SMF_sample, fill_value='extrapolate')(x[k])
        sigma2 = yerr[k]**2
        result += -0.5 * np.sum((y[k] - y_th) ** 2 / sigma2 + np.log(sigma2))

    return result

def log_likelihood_interpolated(x, y, yerr):
    log_par1_span = np.linspace(2,8,150)
    
    result = np.array(progress_map(SMF_single, log_par1_span, n_cpu=None))
   
    interpolated_likelihood = scipy.interpolate.interp1d(log_par1_span, result, fill_value='extrapolate')
    return interpolated_likelihood

#log_likelihood_int = log_likelihood_interpolated(x, y, yerr)
#with open('Puebla_SMF_nDGP_likelihood.pkl', 'wb') as f:
#    pickle.dump(log_likelihood_int, f)
with open('Puebla_SMF_nDGP_likelihood.pkl', 'rb') as f:
    log_likelihood_int = pickle.load(f)

def log_likelihood(theta, x, y, yerr):
    log_par1 = theta
    return log_likelihood_int(log_par1)
    
def log_prior(theta):
    log_par1 = theta
    if (2 < log_par1 < 8):
        return 0
    return -np.inf

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)



nwalkers = 50
ndim = 1
pool_cpu = Pool(8)
sampler = emcee.EnsembleSampler(
    nwalkers, ndim, log_probability, args=(x, y, yerr), pool = pool_cpu
)

initial_params = [3]
per = 0.5
initial_pos = [initial_params + per * np.random.randn(ndim) for _ in range(nwalkers)]
sampler.run_mcmc(initial_pos, 100, progress=True)

flat_samples = sampler.get_chain(flat=True)

import getdist
from getdist import plots, MCSamples

labels = [r'$\log_{10}r_c$']
names = [r'$\log_{10}r_c$']
samples = MCSamples(samples=flat_samples,names = names, labels = labels)
samples.saveAsText("Puebla_nDGP_SMF")
# 1D marginalized comparison plot
g = plots.get_subplot_plotter()
g.plot_1d([samples], r'$\log_{10}r_c$')



plt.savefig('mcmc.pdf')