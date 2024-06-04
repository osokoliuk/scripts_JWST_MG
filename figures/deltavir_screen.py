import sys
sys.path.insert(0, "../")

from JWST_MG.UVLF import UVLF
from JWST_MG.HMF import HMF
from JWST_MG.reionization import reionization
from JWST_MG.delta_c import delta_c
from multiprocessing import Pool, Queue

from JWST_MG.constants import *
from multiprocessing import Queue, Process, cpu_count
import numpy as np
import matplotlib.pyplot as plt
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
from JWST_MG.UVLF import UVLF
from JWST_MG.HMF import HMF
from JWST_MG.reionization import reionization
from JWST_MG.delta_c import delta_c

from JWST_MG.constants import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.rcParams.update({"text.usetex": True})


class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(
            self, useOffset=offset, useMathText=mathText)

    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom

    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
            self.format = r'$\mathdefault{%s}$' % self.format


plt.cla()
plt.figure()
plt.rcParams.update({"text.usetex": True})
fig = plt.figure(figsize=(4.25*1*.95, 2*2*1.05))


ax = plt.subplot(2, 1, 1)
ax.xaxis.set_ticks([0, 0.25, 0.5, 0.75, 1])


ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())


plt.tick_params(axis='both', which='major', direction="in",
                labelsize=14, length=5, top=True, right=True)
plt.tick_params(axis='both', which='minor', direction="in",
                labelsize=11, length=4, top=True, right=True)
plt.tick_params(axis='both', which='major',
                direction="in", labelsize=14, length=5)
plt.tick_params(axis='both', which='minor',
                direction="in", labelsize=11, length=4)


model = 'nDGP'
model_H = 'nDGP'
model_SFR = 'toy'

pars1 = np.logspace(3, 5, 10)
ac_arr = np.linspace(0.05, 1, 20)
par2 = 0

n = len(pars1)
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#8da0cb"])

colors = cmap3(np.linspace(0, 1, n))

import matplotlib.colors as colorss

for i in range(len(pars1)):
    par1 = pars1[i]
    reion_arr = [reionization(np.linspace(ai, ac, 10000), model, model_H, par1, par2) for ac in ac_arr]
    
    def reion(i, model,model_H,par1,par2,a_arr):
        return reion_arr[i].Delta_vir(model,model_H,par1,par2,a_arr)
    pool_cpu = Pool(8)
    iterable = [(i, model,model_H,par1,par2,np.linspace(ai, ac, 10000)) for i, ac in enumerate(ac_arr)]
    a_vir, Deltavir, a_arr, mu_arr = zip(*pool_cpu.starmap(reion,tqdm(iterable, total=len(ac_arr))))
    print(Deltavir)
    plt.plot(ac_arr, Deltavir, c=colors[i], lw=1)


norm = colorss.LogNorm(pars1.min(), pars1.max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap3, norm=norm), ax=ax)
cbar.set_label(r'$r_c$', fontsize=16)

plt.ylabel(r'$\Delta_{\rm vir}(a_c)$', size='16')

# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")
ax.set_xticklabels([])

plt.axhline(18*np.pi**2, c='tab:gray', lw=0.8)
plt.text(0.72, 191, r'$\Delta_{\rm vir}|_{G_{\rm eff}=1}$',
         fontsize=11, c='tab:grey')
ax.fill_between([-0.1, 1.1], 18*np.pi**2-2.5, 18*np.pi**2 +
                2.5, alpha=0.25, color='tab:gray')
plt.xlim(0.0, 1.1)


ax = plt.subplot(2, 1, 2)
ax.xaxis.set_ticks([0, 0.25, 0.5, 0.75, 1])


ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())


plt.tick_params(axis='both', which='major', direction="in",
                labelsize=14, length=5, top=True, right=True)
plt.tick_params(axis='both', which='minor', direction="in",
                labelsize=11, length=4, top=True, right=True)
plt.tick_params(axis='both', which='major',
                direction="in", labelsize=14, length=5)
plt.tick_params(axis='both', which='minor',
                direction="in", labelsize=11, length=4)


model = 'kmoufl'
model_H = 'kmoufl'
model_SFR = 'toy'

pars2 = np.linspace(0.0, 1, 10)
ac_arr = np.linspace(0.05, 1, 20)
pars1 = np.array([0.1, 0.3, 0.5])

n = len(pars2)
cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#66c2a5"]) 
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#fc8d62"]) 
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#8da0cb"])

colors = np.array([cmap1(np.linspace(0, 1, n)), cmap2(np.linspace(0, 1, n)), cmap3(np.linspace(0, 1, n))])

for j in range(len(pars1)):
    par1 = pars1[j]
    for i in range(len(pars2)):
        par2 = pars2[i]
        reion_arr = [reionization(np.linspace(ai, ac, 10000), model, model_H, par1, par2) for ac in ac_arr]
        
        def reion(i, model,model_H,par1,par2,a_arr):
            return reion_arr[i].Delta_vir(model,model_H,par1,par2,a_arr)
        pool_cpu = Pool(8)
        iterable = [(i, model,model_H,par1,par2,np.linspace(ai, ac, 10000)) for i, ac in enumerate(ac_arr)]
        a_vir, Deltavir = zip(*pool_cpu.starmap(reion,tqdm(iterable, total=len(ac_arr))))
    
        plt.plot(ac_arr, Deltavir, c=colors[j][i], lw=1, alpha =0.5)


norm = plt.Normalize(pars2.min(), pars2.max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Grays, norm=norm), ax=ax)
cbar.set_label(r'$K_0$', fontsize=16)

plt.ylabel(r'$\Delta_{\rm vir}(a_c)$', size='16')
plt.xlabel(r'$a_c$', size=16)
# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")

h, l = ax.get_legend_handles_labels()


line1 = Line2D([0], [0], label=r'$\beta=0.1$', color='#66c2a5')
line2 = Line2D([0], [0], label=r'$\beta=0.3$', color='#fc8d62')
line3 = Line2D([0], [0], label=r'$\beta=0.5$', color='#8da0cb')
h.extend([line1, line2, line3])
kw = dict(ncol=1,
          fancybox=True, fontsize=10, frameon=False)
# leg1 = ax.legend(h[:], l[:], bbox_to_anchor=[0.5, 1.08], **kw)
ax.legend(handles=h, loc='upper left', **kw)
plt.axhline(18*np.pi**2, c='tab:gray', lw=0.8)
plt.text(0.72, 191, r'$\Delta_{\rm vir}|_{G_{\rm eff}=1}$',
         fontsize=11, c='tab:grey')
ax.fill_between([-0.1, 1.1], 18*np.pi**2-2.5, 18*np.pi**2 +
                2.5, alpha=0.25, color='tab:gray')
plt.xlim(0.0, 1.1)
plt.tight_layout()

plt.savefig('Delta_vir_screen.pdf', bbox_inches='tight')
