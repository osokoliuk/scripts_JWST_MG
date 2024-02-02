
from JWST_MG.HMF import HMF
from JWST_MG.constants import *
model = 'nDGP'
model_H = 'nDGP'
model_SFR = 'phenomenological_regular'

par1 = 150000
par2 = 1

Masses = np.logspace(7,16,50)
f0 = 0.05

z_arr = np.array([0])
for i in range(len(z_arr)):
    z = z_arr[i]
    a = 1/(1+z)
    HMF_library = HMF(a, model, model_H, par1, par2, Masses/h)
    HMF_fid = HMF_library.ST_mass_function(rhom, Masses/h, a, model_H, model, par1, par2)
    plt.loglog(Masses, HMF_fid)


from hmf import MassFunction     # The main hmf class
import numpy
from numpy import log, log10, sqrt, exp, pi, arange, logspace, linspace, r_, c_
# Don't judge me: 
from matplotlib.pyplot import *
from scipy.interpolate import UnivariateSpline as US
import hmf
from astropy.cosmology import FlatLambdaCDM

h0=0.6732
om0=0.3158
ob0h2=0.022383
ob0=ob0h2/h0**2
ns=0.96605
sig8=0.8120
t_cmb=2.7255
# f_bary=Omega_b/Omega_m
fbary=ob0h2/h0**2/om0


# set up cosmological model:
planck2020_model=FlatLambdaCDM(H0 = 100*h0, Om0=om0, Tcmb0 = t_cmb, Ob0 = ob0)
mf = MassFunction(Mmin= 7, Mmax = 16, hmf_model='SMT', z = 0,
                       cosmo_model=planck2020_model,
                       sigma_8=sig8, n=ns,
                       transfer_model=hmf.density_field.transfer_models.CAMB,
                       transfer_params={'extrapolate_with_eh':True})

Masses= mf.m/h0
HMF = mf.dndm*h0**4
plt.loglog(Masses, HMF)

#my_xlims=np.r_[1e8, 1e12]
#my_ylims=np.r_[1e2, 1e8]
#plt.xlim(*my_xlims)
#plt.ylim(*my_ylims)

plt.savefig('HMF.pdf')
#delta_c_at_ac(self, ac, model, model_H, par1, par2):
