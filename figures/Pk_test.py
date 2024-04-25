# some_file.py
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import odeint
from scipy import integrate
from scipy.special import lambertw
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.pylab as pl
import matplotlib as mpl
import matplotlib
from scipy import interpolate
import scipy
import sys
sys.path.insert(0, "../")
sys.path.insert(0,'../observational_data/GSMF')
from JWST_MG.reionization import reionization
from JWST_MG.cosmological_functions import cosmological_functions
from JWST_MG.constants import *
from JWST_MG.HMF import HMF
from JWST_MG.SMF import SMF

plt.rcParams.update({"text.usetex": True})


model = 'E11'
model_H = 'LCDM'
model_SFR = 'double_power'
pars1 = [1,0]
par2 = 0.0
f0 = 0.12
for i, par1 in enumerate(pars1):
    z_int = z_i = 0
    SMF_library = SMF(1/(1+z_int), model, model_H, model_SFR, par1, par2, 1e8, f0)
    HMF_library = HMF(1/(1+z_i), model, model_H, par1, par2, 1e8)
    k = kvec/h
    Pk = np.array(HMF_library.Pk(1/(1+z_i), model, par1, par2))*h**3
    plt.loglog(k,Pk)

plt.xlim(1e-4,1e0)
plt.ylim(1e2,1e5)
plt.savefig('HMF.pdf')
