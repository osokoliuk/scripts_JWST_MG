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


x = [None]*4
y = [None]*4
yerr = [None]*4

Duncan_z4 = np.loadtxt("/home/oleksii/Downloads/Duncan_z4.dat")
#@plt.errorbar(10**Duncan_z4[:,0],Duncan_z4[:,1],yerr=(Duncan_z4[:,2],Duncan_z4[:,3]), ls = 'None', marker = 's', markercolor = 'tab:orange')
x[0] = 10**Duncan_z4[:,0]
y[0] = Duncan_z4[:,1]
yerr[0] = np.array(Duncan_z4[:,2] + Duncan_z4[:,3])

Duncan_z4 = np.loadtxt("/home/oleksii/Downloads/Duncan_z5.dat")
#@plt.errorbar(10**Duncan_z4[:,0],Duncan_z4[:,1],yerr=(Duncan_z4[:,2],Duncan_z4[:,3]), ls = 'None', marker = 's', markercolor = 'tab:orange')
x[1] = 10**Duncan_z4[:,0]
y[1] = Duncan_z4[:,1]
yerr[1] = np.array(Duncan_z4[:,2] + Duncan_z4[:,3])

Duncan_z4 = np.loadtxt("/home/oleksii/Downloads/Duncan_z6.dat")
#@plt.errorbar(10**Duncan_z4[:,0],Duncan_z4[:,1],yerr=(Duncan_z4[:,2],Duncan_z4[:,3]), ls = 'None', marker = 's', markercolor = 'tab:orange')
x[2] = 10**Duncan_z4[:,0]
y[2] = Duncan_z4[:,1]
yerr[2] = np.array(Duncan_z4[:,2] + Duncan_z4[:,3])

Duncan_z4 = np.loadtxt("/home/oleksii/Downloads/Duncan_z7.dat")
#@plt.errorbar(10**Duncan_z4[:,0],Duncan_z4[:,1],yerr=(Duncan_z4[:,2],Duncan_z4[:,3]), ls = 'None', marker = 's', markercolor = 'tab:orange')
x[3] = 10**Duncan_z4[:,0]
y[3] = Duncan_z4[:,1]
yerr[3] = np.array(Duncan_z4[:,2] + Duncan_z4[:,3])

model = 'nDGP'
model_H = 'nDGP'
model_SFR = 'Behroozi'
par2 = 1
f0 = 0.3

def SMF_func(z, par1):
    SMF_library = SMF(1/(1+z), model, model_H, model_SFR, par1, par2, 1e8, f0)
    HMF_library = HMF(1/(1+z), model, model_H, par1, par2, 1e8)
    Pk = np.array(HMF_library.Pk(1/(1+z), model, par1, par2))*h**3
    k = kvec/h
    Masses = np.logspace(8,14,100)
    Masses_star, SMF_sample = SMF_library.SMF_obs(Masses, rhom, 1/(1+z), model_H, model, model_SFR, par1, par2, k, Pk, f0)
    return Masses_star, SMF_sample

def log_likelihood(theta, x, y, yerr):
    log_par1 = theta
    if hasattr(log_par1, '__len__') and (not isinstance(log_par1, str)):
        result = np.zeros(len(log_par1), dtype=np.float64)
        for i, log_par1 in enumerate(log_par1):
            par1 = 10**log_par1
            Masses_star, SMF_sample = SMF_func(4, par1)
            y_th = scipy.interpolate.interp1d(Masses_star, SMF_sample, fill_value='extrapolate')(x[0])
            sigma2 = yerr[0]**2
            result[i] += -0.5 * np.sum((y[0] - y_th) ** 2 / sigma2 + np.log(sigma2))

            Masses_star, SMF_sample = SMF_func(5, par1)
            y_th = scipy.interpolate.interp1d(Masses_star, SMF_sample, fill_value='extrapolate')(x[1])
            sigma2 = yerr[1]**2
            result[i] += -0.5 * np.sum((y[1] - y_th) ** 2 / sigma2 + np.log(sigma2))

            Masses_star, SMF_sample = SMF_func(6, par1)
            y_th = scipy.interpolate.interp1d(Masses_star, SMF_sample, fill_value='extrapolate')(x[2])
            sigma2 = yerr[2]**2
            result[i] += -0.5 * np.sum((y[2] - y_th) ** 2 / sigma2 + np.log(sigma2))

            Masses_star, SMF_sample = SMF_func(7, par1)
            y_th = scipy.interpolate.interp1d(Masses_star, SMF_sample, fill_value='extrapolate')(x[3])
            sigma2 = yerr[3]**2
            result[i] += -0.5 * np.sum((y[3] - y_th) ** 2 / sigma2 + np.log(sigma2))

        return np.array(result)

    
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


import emcee

nwalkers = 4
ndim = 1
from multiprocessing import Pool

pool_cpu = Pool(8)

sampler = emcee.EnsembleSampler(
    nwalkers, ndim, log_probability, args=(x, y, yerr), pool = pool_cpu
)

initial_params = [3]
per = 0.1
initial_pos = [initial_params + per * np.random.randn(ndim) for _ in range(nwalkers)]
sampler.run_mcmc(initial_pos, 50, progress=True)

flat_samples = sampler.get_chain(discard=5, thin=5, flat=True)

import getdist
from getdist import plots, MCSamples

labels = ["log(r_c)"]
names = ["log(r_c)"]
samples = MCSamples(samples=flat_samples,names = names, labels = labels)


model_SFR = 'double_power'
def SMF_func(z, par1):
    SMF_library = SMF(1/(1+z), model, model_H, model_SFR, par1, par2, 1e8, f0)
    HMF_library = HMF(1/(1+z), model, model_H, par1, par2, 1e8)
    Pk = np.array(HMF_library.Pk(1/(1+z), model, par1, par2))*h**3
    k = kvec/h
    Masses = np.logspace(8,14,100)
    Masses_star, SMF_sample = SMF_library.SMF_obs(Masses, rhom, 1/(1+z), model_H, model, model_SFR, par1, par2, k, Pk, f0)
    return Masses_star, SMF_sample

def log_likelihood(theta, x, y, yerr):
    log_par1 = theta
    if hasattr(log_par1, '__len__') and (not isinstance(log_par1, str)):
        result = np.zeros(len(log_par1), dtype=np.float64)
        for i, log_par1 in enumerate(log_par1):
            par1 = 10**log_par1
            Masses_star, SMF_sample = SMF_func(4, par1)
            y_th = scipy.interpolate.interp1d(Masses_star, SMF_sample, fill_value='extrapolate')(x[0])
            sigma2 = yerr[0]**2
            result[i] += -0.5 * np.sum((y[0] - y_th) ** 2 / sigma2 + np.log(sigma2))

            Masses_star, SMF_sample = SMF_func(5, par1)
            y_th = scipy.interpolate.interp1d(Masses_star, SMF_sample, fill_value='extrapolate')(x[1])
            sigma2 = yerr[1]**2
            result[i] += -0.5 * np.sum((y[1] - y_th) ** 2 / sigma2 + np.log(sigma2))

            Masses_star, SMF_sample = SMF_func(6, par1)
            y_th = scipy.interpolate.interp1d(Masses_star, SMF_sample, fill_value='extrapolate')(x[2])
            sigma2 = yerr[2]**2
            result[i] += -0.5 * np.sum((y[2] - y_th) ** 2 / sigma2 + np.log(sigma2))

            Masses_star, SMF_sample = SMF_func(7, par1)
            y_th = scipy.interpolate.interp1d(Masses_star, SMF_sample, fill_value='extrapolate')(x[3])
            sigma2 = yerr[3]**2
            result[i] += -0.5 * np.sum((y[3] - y_th) ** 2 / sigma2 + np.log(sigma2))

        return np.array(result)

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


import emcee

nwalkers = 4
ndim = 1
from multiprocessing import Pool

pool_cpu = Pool(8)

sampler = emcee.EnsembleSampler(
    nwalkers, ndim, log_probability, args=(x, y, yerr), pool = pool_cpu
)

initial_params = [3]
per = 0.1
initial_pos = [initial_params + per * np.random.randn(ndim) for _ in range(nwalkers)]
sampler.run_mcmc(initial_pos, 50, progress=True)

flat_samples = sampler.get_chain(discard=5, thin=5, flat=True)

import getdist
from getdist import plots, MCSamples

labels = ["log(r_c)"]
names = ["log(r_c)"]
samples2 = MCSamples(samples=flat_samples,names = names, labels = labels)

# 1D marginalized comparison plot
g = plots.get_single_plotter(width_inch=3)
g.plot_1d([samples, samples2], r'$\log_{10}r_c$')



plt.savefig('mcmc.pdf')