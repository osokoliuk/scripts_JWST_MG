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
fig = plt.figure(figsize=(4.25*1*.95, 4*2*1.25))


ax = plt.subplot(4, 1, 1)

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


model = 'E11'
model_H = 'LCDM'
model_SFR = 'toy'

pars1 = np.linspace(-1, 1, 10)
ac_arr = np.linspace(0.15, 1, 20)
par2 = 0

n = len(pars1)
colors = pl.cm.Blues(np.linspace(0, 1, n))

for i in tqdm(range(len(pars1))):
    Delta = []
    par1 = pars1[i]
    for ac in ac_arr:
        deltac = delta_c(ac, model, model_H, par1, par2)
        dc = deltac.delta_c_at_ac(ac, model, model_H, par1, par2)
        Delta.append(dc)
        # print(Deltavir)
    plt.plot(ac_arr, Delta, c=colors[i], lw=1)


norm = plt.Normalize(pars1.min(), pars1.max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Blues, norm=norm), ax=ax)
cbar.set_label(r'$E_{11}$', fontsize=16)

plt.ylabel(r'$\delta_{\rm l}(a_c)$', size='16')

# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")
ax.set_xticklabels([])


plt.xlim(0.1, 1.1)


ax = plt.subplot(4, 1, 2)

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


model = 'gmu'
model_H = 'LCDM'
model_SFR = 'toy'

pars1 = np.linspace(0, 1, 10)
ac_arr = np.linspace(0.15, 1, 20)
par2 = 0

n = len(pars1)
colors = pl.cm.Blues(np.linspace(0, 1, n))

for i in tqdm(range(len(pars1))):
    Delta = []
    par1 = pars1[i]
    for ac in ac_arr:
        deltac = delta_c(ac, model, model_H, par1, par2)
        dc = deltac.delta_c_at_ac(ac, model, model_H, par1, par2)
        Delta.append(dc)
        # print(Deltavir)
    plt.plot(ac_arr, Delta, c=colors[i], lw=1)


norm = plt.Normalize(pars1.min(), pars1.max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Blues, norm=norm), ax=ax)
cbar.set_label(r'$g_\mu$', fontsize=16)

plt.ylabel(r'$\delta_{\rm l}(a_c)$', size='16')

# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")
ax.set_xticklabels([])


plt.xlim(0.1, 1.1)


ax = plt.subplot(4, 1, 3)

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


model = 'DES'
model_H = 'LCDM'
model_SFR = 'toy'

pars2 = np.linspace(0, 1, 10)
ac_arr = np.linspace(0.15, 1, 20)
pars1 = np.array([0, 1, -1])

n = len(pars2)
colors = np.array([pl.cm.Blues(np.linspace(0, 1, n)), pl.cm.Reds(
    np.linspace(0, 1, n)), pl.cm.Purples(np.linspace(0, 1, n))])

for j in tqdm(range(len(pars1))):
    par1 = pars1[j]
    for i in range(len(pars2)):
        Delta = []
        par2 = pars2[i]
        for ac in ac_arr:
            deltac = delta_c(ac, model, model_H, par1, par2)
            dc = deltac.delta_c_at_ac(ac, model, model_H, par1, par2)
            Delta.append(dc)
            # print(Deltavir)
        plt.plot(ac_arr, Delta, c=colors[j][i], alpha=0.5, lw=1)


norm = plt.Normalize(pars2.min(), pars2.max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Grays, norm=norm), ax=ax)
cbar.set_label(r'$T_1$', fontsize=16)

plt.ylabel(r'$\delta_{\rm l}(a_c)$', size='16')

# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")
ax.set_xticklabels([])

h, l = ax.get_legend_handles_labels()


line1 = Line2D([0], [0], label=r'$T_2=0$', color='tab:blue')
line2 = Line2D([0], [0], label=r'$T_2=1$', color='tab:red')
line3 = Line2D([0], [0], label=r'$T_2=-1$', color='tab:purple')
h.extend([line1, line2, line3])
kw = dict(ncol=1,
          fancybox=True, fontsize=10, frameon=False)
# leg1 = ax.legend(h[:], l[:], bbox_to_anchor=[0.5, 1.08], **kw)
ax.legend(handles=h, **kw)

plt.xlim(0.1, 1.1)


ax = plt.subplot(4, 1, 4)

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


model = 'wCDM'
model_H = 'wCDM'
model_SFR = 'toy'

pars2 = np.linspace(0.4, 0.6, 10)
ac_arr = np.linspace(0.15, 1, 20)
pars1 = np.array([-1.5, -1, -0.5])

n = len(pars2)
colors = np.array([pl.cm.Blues(np.linspace(0, 1, n)), pl.cm.Reds(
    np.linspace(0, 1, n)), pl.cm.Purples(np.linspace(0, 1, n))])

for j in tqdm(range(len(pars1))):
    par1 = pars1[j]
    for i in range(len(pars2)):
        Delta = []
        par2 = pars2[i]
        for ac in ac_arr:
            deltac = delta_c(ac, model, model_H, par1, par2)
            dc = deltac.delta_c_at_ac(ac, model, model_H, par1, par2)
            Delta.append(dc)
            # print(Deltavir)
        plt.plot(ac_arr, Delta, c=colors[j][i], alpha=0.5, lw=1)


norm = plt.Normalize(pars2.min(), pars2.max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Grays, norm=norm), ax=ax)
cbar.set_label(r'$\gamma$', fontsize=16)

plt.ylabel(r'$\delta_{\rm l}(a_c)$', size='16')
plt.xlabel(r'$a_c$', size=16)
# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")

h, l = ax.get_legend_handles_labels()


line1 = Line2D([0], [0], label=r'$w_\Lambda=-1.5$', color='tab:blue')
line2 = Line2D([0], [0], label=r'$w_\Lambda=-1$', color='tab:red')
line3 = Line2D([0], [0], label=r'$w_\Lambda=-0.5$', color='tab:purple')
h.extend([line1, line2, line3])
kw = dict(ncol=1,
          fancybox=True, fontsize=10, frameon=False)
# leg1 = ax.legend(h[:], l[:], bbox_to_anchor=[0.5, 1.08], **kw)
ax.legend(handles=h, **kw)

plt.xlim(0.1, 1.1)


"""
par2 = 0.5
ac_arr = np.array([0.1, 0.5, 0.75, 1])
pars1 = np.linspace(-1.5, -0.5, 15)

n = len(ac_arr)
colors = pl.cm.Reds(np.linspace(0, 1, n))


for i in range(len(ac_arr)):
    Delta = []
    ac = ac_arr[i]
    a_arr = np.linspace(ai, ac, 10000)
    for par1 in pars1:
        reion = reionization(a_arr, model, model_H, par1, par2)
        a_vir, Deltavir = reion.Delta_vir(model, model_H, par1, par2, a_arr)
        print(Deltavir)
        Delta.append(Deltavir)
        # print(Deltavir)
    plt.plot(pars1, Delta, c=color[i][0], ls=line[i])


par2 = 0.4
ac_arr = np.array([0.1, 0.5, 0.75, 1])
pars1 = np.linspace(-1.5, -0.5, 15)
colors = pl.cm.Purples(np.linspace(0, 1, n))



for i in range(len(ac_arr)):
    Delta = []
    ac = ac_arr[i]
    a_arr = np.linspace(ai, ac, 10000)
    for par1 in pars1:
        reion = reionization(a_arr, model, model_H, par1, par2)
        a_vir, Deltavir = reion.Delta_vir(model, model_H, par1, par2, a_arr)
        print(Deltavir)
        Delta.append(Deltavir)
        # print(Deltavir)
    plt.plot(pars1, Delta, c=color[i][0], ls=line[i])
"""

"""
ax = plt.subplot(3, 3, 3)


model = 'DES'
model_H = 'LCDM'
model_SFR = 'toy'

ac_arr = np.linspace(0.1, 1, 10)
pars1 = np.linspace(0, 1, 25)
par2 = 0

n = len(ac_arr)
colors = np.array([pl.cm.Blues(np.linspace(0, 1, n)),  pl.cm.Reds(
    np.linspace(0, 1, n)),  pl.cm.Purples(np.linspace(0, 1, n))])

for i in range(len(ac_arr)):
    ac = ac_arr[i]
    a_arr = np.linspace(ai, ac, 10000)
    Delta = []
    for par1 in pars1:
        reion = reionization(a_arr, model, model_H, par1, par2)
        a_vir, Deltavir = reion.Delta_vir(
            model, model_H, par1, par2, a_arr)
        print(Deltavir)
        Delta.append(Deltavir)
        # print(Deltavir)
    plt.plot(pars1, Delta, c=colors[i][0])

plt.title(r"$T_2=0$")
norm = plt.Normalize(ac_arr.min(), ac_arr.max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Blues, norm=norm), ax=ax)
cbar.set_label(r'$a_{\rm c}$', fontsize=16)

plt.ylabel(r'$\Delta_{\rm vir}(a_{\rm c})$', size='16')
plt.xlabel(r'$T_{1}$', size='16')

# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")

h, l = ax.get_legend_handles_labels()
line1 = Line2D([0], [0], label=r'$T_2=0$', color='tab:blue')
line2 = Line2D([0], [0], label=r'$T_2=5$', color='tab:red')
line3 = Line2D([0], [0], label=r'$T_2=-5$', color='tab:purple')
h.extend([line1, line2, line3])
kw = dict(ncol=3, loc="lower center",
          fancybox=True, fontsize=11, frameon=False)
leg1 = ax.legend(h[:], l[:], bbox_to_anchor=[0.5, 1.08], **kw)
ax.add_artist(leg1)

plt.tight_layout()
"""

plt.savefig('deltac_pheno.pdf', bbox_inches='tight')
