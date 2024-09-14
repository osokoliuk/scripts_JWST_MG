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
mpl.rcParams['axes.linewidth'] = 1.5


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
                direction="in", labelsize=14, length=5, width = 1.5)
plt.tick_params(axis='both', which='minor',
                direction="in", labelsize=11, length=4, width = 1.1)


model = 'E11'
model_H = 'LCDM'
model_SFR = 'toy'

pars1 = np.logspace(-1,1, 10)
ac_arr = np.linspace(0.05, 1, 20)
par2 = 0
n = len(pars1)
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#48639e"])

colors = cmap3(np.linspace(0, 1, n))
for i in tqdm(range(len(pars1))):
    par1 = pars1[i]
    reion_arr = [reionization(np.linspace(ai, ac, 10000), model, model_H, par1, par2) for ac in ac_arr]
    
    def reion(i, model,model_H,par1,par2,a_arr):
        return reion_arr[i].Delta_vir(model,model_H,par1,par2,a_arr)
    pool_cpu = Pool(8)
    iterable = [(i, model,model_H,par1,par2,np.linspace(ai, ac, 10000)) for i, ac in enumerate(ac_arr)]
    a_vir, Deltavir = zip(*pool_cpu.starmap(reion,tqdm(iterable, total=len(ac_arr))))
    print(Deltavir)
    plt.plot(ac_arr, Deltavir, c=colors[i], lw=1)


norm = plt.Normalize(np.log10(pars1).min(), np.log10(pars1).max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap3, norm=norm), ax=ax, pad = 0)
cbar.ax.tick_params(width=1.5, length=5, which = 'major')
cbar.ax.tick_params(width=1.1, length=4, which = 'minor')
cbar.set_label(r'$E_{11}$', fontsize=16)
labels = cbar.ax.get_yticklabels()
cbar.ax.set_yticklabels(labels)
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
# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")


plt.xlim(0, 1.1)


ax = plt.subplot(4, 1, 2)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())


plt.tick_params(axis='both', which='major', direction="in",
                labelsize=14, length=5, top=True, right=True)
plt.tick_params(axis='both', which='minor', direction="in",
                labelsize=11, length=4, top=True, right=True)
plt.tick_params(axis='both', which='major',
                direction="in", labelsize=14, length=5, width = 1.5)
plt.tick_params(axis='both', which='minor',
                direction="in", labelsize=11, length=4, width = 1.1)


model = 'gmu'
model_H = 'LCDM'
model_SFR = 'toy'

pars1 = np.logspace(-1,1, 10)
ac_arr = np.linspace(0.05, 1, 20)
par2 = 0
n = len(pars1)
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#48639e"])

colors = cmap3(np.linspace(0, 1, n))
for i in tqdm(range(len(pars1))):
    par1 = pars1[i]
    reion_arr = [reionization(np.linspace(ai, ac, 10000), model, model_H, par1, par2) for ac in ac_arr]
    
    def reion(i, model,model_H,par1,par2,a_arr):
        return reion_arr[i].Delta_vir(model,model_H,par1,par2,a_arr)
    pool_cpu = Pool(8)
    iterable = [(i, model,model_H,par1,par2,np.linspace(ai, ac, 10000)) for i, ac in enumerate(ac_arr)]
    a_vir, Deltavir = zip(*pool_cpu.starmap(reion,tqdm(iterable, total=len(ac_arr))))
    print(Deltavir)
    plt.plot(ac_arr, Deltavir, c=colors[i], lw=1)


norm = plt.Normalize(np.log10(pars1).min(), np.log10(pars1).max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap3, norm=norm), ax=ax, pad = 0)
cbar.ax.tick_params(width=1.5, length=5, which = 'major')
cbar.ax.tick_params(width=1.1, length=4, which = 'minor')
cbar.set_label(r'$g_{\mu}$', fontsize=16)
labels = cbar.ax.get_yticklabels()
labels[-1] = ""
cbar.ax.set_yticklabels(labels)
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
# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")


plt.xlim(0, 1.1)




ax = plt.subplot(4, 1, 3)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())


plt.tick_params(axis='both', which='major', direction="in",
                labelsize=14, length=5, top=True, right=True)
plt.tick_params(axis='both', which='minor', direction="in",
                labelsize=11, length=4, top=True, right=True)
plt.tick_params(axis='both', which='major',
                direction="in", labelsize=14, length=5, width = 1.5)
plt.tick_params(axis='both', which='minor',
                direction="in", labelsize=11, length=4, width = 1.1)


model = 'DES'
model_H = 'LCDM'
model_SFR = 'toy'

pars2 = np.linspace(-0.75,1,10)
ac_arr = np.linspace(0.05, 1, 10)
pars1 = np.array([-1,0,1])

n = len(pars2)
cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#398e73"]) 
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#e64304"]) 
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#48639e"])

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
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Grays, norm=norm), ax=ax, pad = 0)
cbar.ax.tick_params(width=1.5, length=5, which = 'major')
cbar.ax.tick_params(width=1.1, length=4, which = 'minor')
cbar.set_label(r'$T_1$', fontsize=16)
labels = cbar.ax.get_yticklabels()
labels[-1] = ""
cbar.ax.set_yticklabels(labels)
plt.ylabel(r'$\Delta_{\rm vir}(a_c)$', size='16')
plt.xlabel(r'$a_c$', size=16)
# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")

h, l = ax.get_legend_handles_labels()

plt.axhline(18*np.pi**2, c='tab:gray', lw=0.8)
plt.text(0.72, 191, r'$\Delta_{\rm vir}|_{G_{\rm eff}=1}$',
         fontsize=11, c='tab:grey')
ax.fill_between([-0.1, 1.1], 18*np.pi**2-2.5, 18*np.pi**2 +
                2.5, alpha=0.25, color='tab:gray')





h, l = ax.get_legend_handles_labels()



line1 = Line2D([0], [0], label=r'$T_2=-1$', color='#398e73')
line2 = Line2D([0], [0], label=r'$T_2=0$', color='#e64304')
line3 = Line2D([0], [0], label=r'$T_2=1$', color='#48639e')


h.extend([line1, line2, line3])
kw = dict(ncol=1,
          fancybox=True, fontsize=10, frameon=False)
# leg1 = ax.legend(h[:], l[:], bbox_to_anchor=[0.5, 1.08], **kw)
ax.legend(handles=h, **kw)

plt.xlim(0, 1.1)




ax = plt.subplot(4, 1, 4)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())


plt.tick_params(axis='both', which='major', direction="in",
                labelsize=14, length=5, top=True, right=True)
plt.tick_params(axis='both', which='minor', direction="in",
                labelsize=11, length=4, top=True, right=True)
plt.tick_params(axis='both', which='major',
                direction="in", labelsize=14, length=5, width = 1.5)
plt.tick_params(axis='both', which='minor',
                direction="in", labelsize=11, length=4, width = 1.1)


model = 'wCDM'
model_H = 'wCDM'
model_SFR = 'toy'

pars2 = np.linspace(0.4,0.6,10)
ac_arr = np.linspace(0.05, 1, 10)
pars1 = np.array([-1.5,-1,-0.5])

n = len(pars2)
cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#398e73"]) 
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#e64304"]) 
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#48639e"])

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
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Grays, norm=norm), ax=ax, pad = 0)
cbar.ax.tick_params(width=1.5, length=5, which = 'major')
cbar.ax.tick_params(width=1.1, length=4, which = 'minor')
cbar.set_label(r'$\gamma$', fontsize=16)
labels = cbar.ax.get_yticklabels()
labels[-1] = ""
cbar.ax.set_yticklabels(labels)
plt.ylabel(r'$\Delta_{\rm vir}(a_c)$', size='16')
plt.xlabel(r'$a_c$', size=16)
# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")

h, l = ax.get_legend_handles_labels()

plt.axhline(18*np.pi**2, c='tab:gray', lw=0.8)
plt.text(0.72, 191, r'$\Delta_{\rm vir}|_{G_{\rm eff}=1}$',
         fontsize=11, c='tab:grey')
ax.fill_between([-0.1, 1.1], 18*np.pi**2-2.5, 18*np.pi**2 +
                2.5, alpha=0.25, color='tab:gray')





h, l = ax.get_legend_handles_labels()



line1 = Line2D([0], [0], label=r'$w_\Lambda=-1.5$', color='#398e73')
line2 = Line2D([0], [0], label=r'$w_\Lambda=-1$', color='#e64304')
line3 = Line2D([0], [0], label=r'$w_\Lambda=-0.5$', color='#48639e')


h.extend([line1, line2, line3])
kw = dict(ncol=1,
          fancybox=True, fontsize=10, frameon=False)
# leg1 = ax.legend(h[:], l[:], bbox_to_anchor=[0.5, 1.08], **kw)
ax.legend(handles=h, **kw)

plt.xlim(0, 1.1)



plt.subplots_adjust(wspace=0, hspace=0)

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
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Blues, norm=norm), ax=ax, pad = 0)
cbar.ax.tick_params(width=1.5, length=5, which = 'major')
cbar.ax.tick_params(width=1.1, length=4, which = 'minor')
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

plt.savefig('deltavir_pheno.pdf', bbox_inches='tight')