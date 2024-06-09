import sys
sys.path.insert(0, "../")

from JWST_MG.UVLF import UVLF
from JWST_MG.SMD import SMD
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
import matplotlib.colors as colorss
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
fig = plt.figure(figsize=(4.25*2*.95, 2*2*1.05))


ax = plt.subplot(2, 2, 1)

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
model_SFR = 'Puebla'

pars1 = np.logspace(2.25, 5, 10)
par2 = 0
z = 8
f0 = 0
n = len(pars1)
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#48639e"])

colors = cmap3(np.linspace(0, 1, n))
SMD_library = SMD(1/(1+z), model, model_H, model_SFR, pars1, par2, 1e8, f0)
Pk_arr = []
for par1 in pars1:
    HMF_library = HMF(1/(1+z), model, model_H, par1, par2, 1e8)
    Pk_arr.append(np.array(HMF_library.Pk(1/(1+z), model, par1, par2))*h**3)
k = kvec/h
Masses = np.logspace(8,16,1000)

pool_cpu = Pool(8)
iterable = [(Masses, rhom, 1/(1+z), model_H, model, model_SFR, par1, par2, k, Pk_arr[i], f0) for i,par1 in enumerate(pars1)]
Masses_stars, SMDs = zip(*pool_cpu.starmap(SMD_library.SMD,tqdm(iterable, total=len(pars1))))


#print(SMF)
#Masses_star = SMF[0]
#SMF_obs = SMF[1]
for i, Masses_star in enumerate(Masses_stars):
    plt.loglog(Masses_stars[i], SMDs[i], c = colors[i], lw=1)





data = np.loadtxt('/home/oleksii/codes/scripts_JWST_MG/observational_data/SMD/JWST_z8.txt')
x = data[:,0]
y = data[:,1]
ye_low = data[:,2]
ye_upp =  data[:,3]


plt.errorbar(x,y,yerr=c_[ye_low, ye_upp].T,capsize=0,ecolor='k',color='w',marker='o',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='k')



norm = colorss.LogNorm(pars1.min(), pars1.max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap3, norm=norm), ax=ax)
cbar.set_label(r'$r_c$', fontsize=16)

plt.ylabel(r'$\rho_\star\;[M_\odot\;\rm Mpc^{-3}]$', size='16')

# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")
ax.set_xticklabels([])
plt.xlim(1e5,10**11.5)
plt.ylim(1e3,10**7.5)
plt.title(r"$z=8$", size = 16)


ax = plt.subplot(2, 2, 2)

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
model_SFR = 'Puebla'

pars1 = np.logspace(2.25, 5, 10)
par2 = 0
z = 9
f0 = 0
n = len(pars1)
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#48639e"])

colors = cmap3(np.linspace(0, 1, n))

SMD_library = SMD(1/(1+z), model, model_H, model_SFR, pars1, par2, 1e8, f0)
Pk_arr = []
for par1 in pars1:
    HMF_library = HMF(1/(1+z), model, model_H, par1, par2, 1e8)
    Pk_arr.append(np.array(HMF_library.Pk(1/(1+z), model, par1, par2))*h**3)
k = kvec/h
Masses = np.logspace(8,16,1000)

pool_cpu = Pool(8)
iterable = [(Masses, rhom, 1/(1+z), model_H, model, model_SFR, par1, par2, k, Pk_arr[i], f0) for i,par1 in enumerate(pars1)]
Masses_stars, SMDs = zip(*pool_cpu.starmap(SMD_library.SMD,tqdm(iterable, total=len(pars1))))


#print(SMF)
#Masses_star = SMF[0]
#SMF_obs = SMF[1]
for i, Masses_star in enumerate(Masses_stars):
    plt.loglog(Masses_stars[i], SMDs[i], c = colors[i], lw=1)


data = np.loadtxt('/home/oleksii/codes/scripts_JWST_MG/observational_data/SMD/JWST_z9.txt')
x = data[:,0]
y = data[:,1]
ye_low = data[:,2]
ye_upp =  data[:,3]


plt.errorbar(x,y,yerr=c_[ye_low, ye_upp].T,capsize=0,ecolor='k',color='w',marker='o',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='k')


norm = colorss.LogNorm(pars1.min(), pars1.max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap3, norm=norm), ax=ax)
cbar.set_label(r'$r_c$', fontsize=16)


# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")
ax.set_xticklabels([])
plt.ylabel(r'$\rho_\star\;[M_\odot\;\rm Mpc^{-3}]$', size='16')

plt.xlim(1e5,10**11.5)
plt.ylim(1e3,10**7.5)
plt.title(r"$z=9$", size = 16)


ax = plt.subplot(2, 2, 3)

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
model_SFR = 'Puebla'

pars2 = np.linspace(0.0, 1, 10)
pars1 = np.array([0.1, 0.3, 0.5])
z = 8
f0 = 0
n = len(pars2)
cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#398e73"]) 
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#e64304"]) 
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#48639e"])

colors = np.array([cmap1(np.linspace(0, 1, n)), cmap2(np.linspace(0, 1, n)), cmap3(np.linspace(0, 1, n))])

for j, par1 in enumerate(pars1):
    SMD_library = SMD(1/(1+z), model, model_H, model_SFR, pars1, par2, 1e8, f0)
    Pk_arr = []
    for par2 in pars2:
        HMF_library = HMF(1/(1+z), model, model_H, par1, par2, 1e8)
        Pk_arr.append(np.array(HMF_library.Pk(1/(1+z), model, par1, par2))*h**3)
    k = kvec/h
    Masses = np.logspace(8,16,1000)

    pool_cpu = Pool(8)
    iterable = [(Masses, rhom, 1/(1+z), model_H, model, model_SFR, par1, par2, k, Pk_arr[i], f0) for i,par2 in enumerate(pars2)]
    Masses_stars, SMDs = zip(*pool_cpu.starmap(SMD_library.SMD,tqdm(iterable, total=len(pars2))))


    #print(SMF)
    #Masses_star = SMF[0]
    #SMF_obs = SMF[1]
    for i, Masses_star in enumerate(Masses_stars):
        plt.loglog(Masses_stars[i], SMDs[i],  c=colors[j][i], alpha=0.5)


hhh, llll = ax.get_legend_handles_labels()


line1 = Line2D([0], [0], label=r'$\beta=0.1$', color='#398e73')
line2 = Line2D([0], [0], label=r'$\beta=0.3$', color='#e64304')
line3 = Line2D([0], [0], label=r'$\beta=0.5$', color='#48639e')
hhh.extend([line1, line2, line3])
kw = dict(ncol=1,
          fancybox=True, fontsize=10, frameon=False)
# leg1 = ax.legend(h[:], l[:], bbox_to_anchor=[0.5, 1.08], **kw)
ax.legend(handles=hhh, **kw,loc='lower left')



data = np.loadtxt('/home/oleksii/codes/scripts_JWST_MG/observational_data/SMD/JWST_z8.txt')
x = data[:,0]
y = data[:,1]
ye_low = data[:,2]
ye_upp =  data[:,3]


plt.errorbar(x,y,yerr=c_[ye_low, ye_upp].T,capsize=0,ecolor='k',color='w',marker='o',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='k')


norm = plt.Normalize(pars2.min(), pars2.max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Grays, norm=norm), ax=ax)
cbar.set_label(r'$K_0$', fontsize=16)

plt.ylabel(r'$\rho_\star\;[M_\odot\;\rm Mpc^{-3}]$', size='16')
plt.xlabel(r'$\rm M_\star\;[M_\odot]$', size = '16')

# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")
plt.xlim(1e5,10**11.5)
plt.ylim(10**2.25,10**7.75)



ax = plt.subplot(2, 2, 4)

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
model_SFR = 'Puebla'

pars2 = np.linspace(0.0, 1, 10)
pars1 = np.array([0.1, 0.3, 0.5])
z = 9
f0 = 0
n = len(pars2)
cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#398e73"]) 
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#e64304"]) 
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#48639e"])

colors = np.array([cmap1(np.linspace(0, 1, n)), cmap2(np.linspace(0, 1, n)), cmap3(np.linspace(0, 1, n))])

for j, par1 in enumerate(pars1):
    SMD_library = SMD(1/(1+z), model, model_H, model_SFR, pars1, par2, 1e8, f0)
    Pk_arr = []
    for par2 in pars2:
        HMF_library = HMF(1/(1+z), model, model_H, par1, par2, 1e8)
        Pk_arr.append(np.array(HMF_library.Pk(1/(1+z), model, par1, par2))*h**3)
    k = kvec/h
    Masses = np.logspace(8,16,1000)

    pool_cpu = Pool(8)
    iterable = [(Masses, rhom, 1/(1+z), model_H, model, model_SFR, par1, par2, k, Pk_arr[i], f0) for i,par2 in enumerate(pars2)]
    Masses_stars, SMDs = zip(*pool_cpu.starmap(SMD_library.SMD,tqdm(iterable, total=len(pars2))))


    #print(SMF)
    #Masses_star = SMF[0]
    #SMF_obs = SMF[1]
    for i, Masses_star in enumerate(Masses_stars):
        plt.loglog(Masses_stars[i], SMDs[i],  c=colors[j][i], alpha=0.5)



hhh, lll = ax.get_legend_handles_labels()

line1 = Line2D([0], [0], label=r'$\beta=0.1$', color='#398e73')
line2 = Line2D([0], [0], label=r'$\beta=0.3$', color='#e64304')
line3 = Line2D([0], [0], label=r'$\beta=0.5$', color='#48639e')
hhh.extend([line1, line2, line3])
kw = dict(ncol=1,
          fancybox=True, fontsize=10, frameon=False)
# leg1 = ax.legend(h[:], l[:], bbox_to_anchor=[0.5, 1.08], **kw)
ax.legend(handles=hhh, **kw,loc='lower left')


data = np.loadtxt('/home/oleksii/codes/scripts_JWST_MG/observational_data/SMD/JWST_z9.txt')
x = data[:,0]
y = data[:,1]
ye_low = data[:,2]
ye_upp =  data[:,3]


plt.errorbar(x,y,yerr=c_[ye_low, ye_upp].T,capsize=0,ecolor='k',color='w',marker='o',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='k')


norm = plt.Normalize(pars2.min(), pars2.max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Grays, norm=norm), ax=ax)
cbar.set_label(r'$K_0$', fontsize=16)

plt.ylabel(r'$\rho_\star\;[M_\odot\;\rm Mpc^{-3}]$', size='16')
plt.xlabel(r'$\rm M_\star\;[M_\odot]$', size = '16')
# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")
plt.xlim(1e5,10**11.5)
plt.ylim(10**2.25,10**7.75)





"""
ax = plt.subplot(2, 1, 2)

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
ac_arr = np.linspace(0.15, 1, 20)
pars1 = np.array([0.1, 0.3, 0.5])

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
        plt.plot(ac_arr, Delta, c=colors[j][i], alpha=1, lw=1)


norm = plt.Normalize(pars2.min(), pars2.max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Grays, norm=norm), ax=ax)
cbar.set_label(r'$K_0$', fontsize=16)

plt.ylabel(r'$\delta_{\rm l}(a_c)$', size='16')

# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")

h, l = ax.get_legend_handles_labels()


line1 = Line2D([0], [0], label=r'$\beta=0.1$', color='tab:blue')
line2 = Line2D([0], [0], label=r'$\beta=0.3$', color='tab:red')
line3 = Line2D([0], [0], label=r'$\beta=0.5$', color='tab:purple')
h.extend([line1, line2, line3])
kw = dict(ncol=1,
          fancybox=True, fontsize=10, frameon=False)
# leg1 = ax.legend(h[:], l[:], bbox_to_anchor=[0.5, 1.08], **kw)
ax.legend(handles=h, **kw,loc='upper right')

plt.xlim(0.1, 1.1)
plt.xlabel(r'$a_c$', size=16)

585 Kyiv comet station. Observers A. Kasianchuk, M. Solomakha, A. Simon, M. Bilodid, O. Pastoven, D. Provolovska, A. Baransky, M. Derkach. Measurers A. Baransky, O. Sokoliuk, L. Sotnichenko, V. Vdovenko, A. Dzygunenko. 0.7-m f/4 reflector + CCD. 
"""

plt.tight_layout()

plt.savefig('SMD_screen_Puebla.pdf', bbox_inches='tight')
