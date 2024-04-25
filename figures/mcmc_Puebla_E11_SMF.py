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


model = 'E11'
model_H = 'LCDM'
model_SFR = 'Puebla'
f0 = 0

def SMF_func(z, par1, par2):
    SMF_library = SMF(1/(1+z), model, model_H, model_SFR, par1, par2, 1e8, f0)
    HMF_library = HMF(1/(1+z), model, model_H, par1, par2, 1e8)
    Pk = np.array(HMF_library.Pk(1/(1+z), model, par1, par2))*h**3
    k = kvec/h
    Masses = np.logspace(6,16,100)
    Masses_star, SMF_sample = SMF_library.SMF_obs(Masses, rhom, 1/(1+z), model_H, model, model_SFR, par1, par2, k, Pk, f0)
    return Masses_star, SMF_sample

def SMF_single(par1, par2):
    result = 0
    for k, zi in enumerate(zs):
        Masses_star, SMF_sample = SMF_func(zi, par1, par2)
        y_th = scipy.interpolate.interp1d(Masses_star, SMF_sample, fill_value='extrapolate')(x[k])
        sigma2 = yerr[k]**2
        result += -0.5 * np.sum((y[k] - y_th) ** 2 / sigma2 + np.log(sigma2))

    return result

from scipy.interpolate import LinearNDInterpolator as linterp
from scipy.interpolate import NearestNDInterpolator as nearest

class LinearNDInterpolatorExt(object):
    def __init__(self, points, values):
        self.funcinterp = linterp(points, values)
        self.funcnearest = nearest(points, values)
    
    def __call__(self, *args):
        z = self.funcinterp(*args)
        chk = np.isnan(z)
        if chk.any():
            return np.where(chk, self.funcnearest(*args), z)
        else:
            return z




def log_likelihood_interpolated(x, y, yerr):
    sampler = qmc.LatinHypercube(d=2)
    sample = sampler.random(n=450)
    l_bounds = [-1,-1]
    u_bounds = [1,1]
    sample_scaled = qmc.scale(sample, l_bounds, u_bounds)
    par1_span = sample_scaled[:,0]
    par2_span = sample_scaled[:,1]
    #xx, yy = np.meshgrid(log_par1_span, f0_span)
    iterable = []
    for par1, par2 in zip(par1_span,par2_span):
        iterable.append([par1,par2])
    #result = zip(*pool_cpu.starmap(SMF_single,tqdm(iterable, total=len(log_par1)*len(f0_span))))

    result = np.array(progress_starmap(SMF_single, iterable, n_cpu=None))

    interpolated_likelihood = LinearNDInterpolatorExt(list(zip(par1_span, par2_span)),result)
    return interpolated_likelihood

#log_likelihood_int = log_likelihood_interpolated(x, y, yerr)
#with open('Puebla_SMF_E11_likelihood.pkl', 'wb') as f:
#    pickle.dump(log_likelihood_int, f)
with open('Puebla_SMF_E11_likelihood.pkl', 'rb') as f:
    log_likelihood_int = pickle.load(f)

def log_likelihood(theta, x, y, yerr):
    par1, par2 = theta
    return log_likelihood_int(par1, par2)


def log_prior(theta):
    par1, par2 = theta
    if (-1 < par1 < 1 and -1 < par2 < 1):
        return 0
    return -np.inf

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)


import emcee

nwalkers = 50
ndim = 2
from multiprocessing import Pool

pool_cpu = Pool(8)

sampler = emcee.EnsembleSampler(
    nwalkers, ndim, log_probability, args=(x, y, yerr), pool = pool_cpu
)

initial_params = [1, 1]
per = 0.01
initial_pos = [initial_params + per * np.random.randn(ndim) for _ in range(nwalkers)]
sampler.run_mcmc(initial_pos, 5000, progress=True)

flat_samples = sampler.get_chain(discard=5, thin=5, flat=True)

import getdist
from getdist import plots, MCSamples

labels = [r'$E_{11}$', r'$E_{22}$']
names = [r'$E_{11}$', r'$E_{22}$']
samples = MCSamples(samples=flat_samples,names = names, labels = labels)
samples.saveAsText("Puebla_E11_SMF")
# 1D marginalized comparison plot
g = plots.get_subplot_plotter()
g.triangle_plot([samples], filled = True)



plt.savefig('mcmc.pdf')