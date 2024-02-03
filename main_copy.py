
from JWST_MG.HMF import HMF
from JWST_MG.constants import *
model = 'nDGP'
model_H = 'nDGP'
model_SFR = 'phenomenological_regular'

par1 = 150000
par2 = 1

Masses = np.logspace(6,13,50)
f0 = 0.05

z_arr = np.array([9.1])
for i in range(len(z_arr)):
    z = z_arr[i]
    a = 1/(1+z)
    HMF_library = HMF(a, model, model_H, par1, par2, Masses)
    HMF_fid = HMF_library.ST_mass_function(rhom, Masses, a, model_H, model, par1, par2)
    HMF_fid = HMF_fid
    plt.loglog(Masses, HMF_fid)


from hmf import MassFunction     # The main hmf class
import numpy
from numpy import log, log10, sqrt, exp, pi, arange, logspace, linspace, r_, c_
# Don't judge me: 
from matplotlib.pyplot import *
from scipy.interpolate import UnivariateSpline as US
import hmf
from astropy.cosmology import FlatLambdaCDM

h0=h
om0=Omegam0
ob0h2=Omegab0*h**2
ob0=ob0h2/h**2
ns=0.96605
sig8=0.8120
t_cmb=2.7255
# f_bary=Omega_b/Omega_m
fbary=ob0h2/h0**2/om0


# set up cosmological model:
planck2020_model=FlatLambdaCDM(H0 = 100*h0, Om0=om0, Tcmb0 = t_cmb, Ob0 = ob0)
mf = MassFunction(Mmin= 6, Mmax = 13, hmf_model='SMT', z = 9.1,
                       cosmo_model=planck2020_model,
                       sigma_8=sig8, n=ns,
                       transfer_model=hmf.density_field.transfer_models.CAMB,
                       transfer_params={'extrapolate_with_eh':True})

Masses= mf.m/h0
HMF = mf.dndm*h0**4
plt.loglog(Masses, HMF)

Masses = np.logspace(6,13,50)

from colossus.cosmology import cosmology
model_SFR = 'double_power'    
cosmology.setCosmology('planck18');
from colossus.lss import mass_function
mfunc = mass_function.massFunction(Masses, z, mdef = 'fof', model = 'sheth99', q_out = 'dndlnM', sigma_args = {'filt': 'tophat'})

plt.plot(Masses/h,(mfunc/Masses)*h**4)

#my_xlims=np.r_[1e8, 1e12]
#my_ylims=np.r_[1e2, 1e8]
#plt.xlim(*my_xlims)
#plt.ylim(*my_ylims)

plt.savefig('HMF.pdf')
#delta_c_at_ac(self, ac, model, model_H, par1, par2):
