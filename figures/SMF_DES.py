

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
import matplotlib.colors as colorss

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
fig = plt.figure(figsize=(4.25*1*.95, 1*2*0.95))
"""
ax = plt.subplot(4, 1, 1)


ax.xaxis.set_ticks([0, 0.25, 0.5, 0.75, 1])
ax.yaxis.set_ticks([0, 0.5, 1, 1.5])


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

par2 = 6/11
ac = 1
a_arr = np.linspace(ai, ac, 10000)
pars1 = np.linspace(1, -1, 2)

n = len(pars1)
colors = pl.cm.Blues(np.linspace(0, 1, n))

for i in range(len(pars1)):
    par1 = pars1[i]
    deltac = delta_c(ac, model, model_H, par1, par2)
    deltai_collapse = deltac.binary_search_di(
        ac, model, model_H, par1, par2, 0, len(delta_ini), abs_err)
    delta = deltac.linear(deltai_collapse, ac,
                          model, model_H, par1, par2)
    print("a")
    # print(Deltavir)
    plt.plot(delta[:, 0], delta[:, 1], c=colors[i])


norm = plt.Normalize(pars1.min(), pars1.max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Blues, norm=norm), ax=ax)
cbar.set_label(r'$T_1$', fontsize=16)

plt.ylabel(r'$\delta_{\rm m}(a)$', size='16')

# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")
ax.set_xticklabels([])

h, l = ax.get_legend_handles_labels()
kw = dict(ncol=3, loc="lower center",
          fancybox=True, fontsize=11, frameon=False)
leg1 = ax.legend(h[:], l[:], bbox_to_anchor=[0.5, 1.08], **kw)
ax.add_artist(leg1)
plt.axhline(1.675, c='tab:gray', lw=0.8)
plt.text(0.1, 1.375, r'$\delta_{\rm c}|_{G_{\rm eff}=1}$',
         fontsize=13, c='tab:grey')
ax.fill_between([-0.1, 1.1], 1.675-0.025, 1.675 +
                0.025, alpha=0.25, color='tab:gray')
plt.xlim(-0.05, 1.05)

ax = plt.subplot(4, 1, 2)


ax.xaxis.set_ticks([0, 0.25, 0.5, 0.75, 1])
ax.yaxis.set_ticks([0, 0.5, 1, 1.5])


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

par2 = 6/11
ac = 1
a_arr = np.linspace(ai, ac, 10000)
pars1 = np.linspace(0, 1, 2)

n = len(pars1)
colors = pl.cm.Blues(np.linspace(0, 1, n))

for i in range(len(pars1)):
    par1 = pars1[i]
    deltac = delta_c(ac, model, model_H, par1, par2)
    deltai_collapse = deltac.binary_search_di(
        ac, model, model_H, par1, par2, 0, len(delta_ini), abs_err)
    delta = deltac.linear(deltai_collapse, ac,
                          model, model_H, par1, par2)
    print("a")
    # print(Deltavir)
    plt.plot(delta[:, 0], delta[:, 1], c=colors[i])


norm = plt.Normalize(pars1.min(), pars1.max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Blues, norm=norm), ax=ax)
cbar.set_label(r'$g_{\mu}$', fontsize=16)

plt.ylabel(r'$\delta_{\rm m}(a)$', size='16')

# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")
ax.set_xticklabels([])

h, l = ax.get_legend_handles_labels()
kw = dict(ncol=3, loc="lower center",
          fancybox=True, fontsize=11, frameon=False)
leg1 = ax.legend(h[:], l[:], bbox_to_anchor=[0.5, 1.08], **kw)
ax.add_artist(leg1)
plt.axhline(1.675, c='tab:gray', lw=0.8)
plt.text(0.1, 1.375, r'$\delta_{\rm c}|_{G_{\rm eff}=1}$',
         fontsize=13, c='tab:grey')
ax.fill_between([-0.1, 1.1], 1.675-0.025, 1.675 +
                0.025, alpha=0.25, color='tab:gray')
plt.xlim(-0.05, 1.05)
"""

"""
ax = plt.subplot(4, 1, 3)


ax.xaxis.set_ticks([0, 0.25, 0.5, 0.75, 1])

ax.yaxis.set_ticks([0, 0.5, 1, 1.5])

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

pars2 = np.array([0, 5, -5])
ac = 1
a_arr = np.linspace(ai, ac, 10000)
pars1 = np.linspace(0, 1, 2)

n = len(pars1)
colors = np.array([pl.cm.Blues(np.linspace(0, 1, n)), pl.cm.Reds(
    np.linspace(0, 1, n)), pl.cm.Purples(np.linspace(0, 1, n))])

for j in range(len(pars2)):
    par2 = pars2[j]
    for i in range(len(pars1)):
        par1 = pars1[i]
        deltac = delta_c(ac, model, model_H, par1, par2)
        deltai_collapse = deltac.binary_search_di(
            ac, model, model_H, par1, par2, 0, len(delta_ini), abs_err)
        delta = deltac.linear(deltai_collapse, ac,
                              model, model_H, par1, par2)
        print("a")
        # print(Deltavir)
        plt.plot(delta[:, 0], delta[:, 1], c=colors[j][i], alpha=0.3, lw=0.75)


norm = plt.Normalize(pars1.min(), pars1.max())
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Blues, norm=norm), ax=ax)
cbar.set_label(r'$T_{1}$', fontsize=16)

plt.ylabel(r'$\delta_{\rm m}(a)$', size='16')

# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")
ax.set_xticklabels([])

h, l = ax.get_legend_handles_labels()


line1 = Line2D([0], [0], label=r'$T_2=0$', color='tab:blue')
line2 = Line2D([0], [0], label=r'$T_2=5$', color='tab:red')
line3 = Line2D([0], [0], label=r'$T_2=-5$', color='tab:purple')
h.extend([line1, line2, line3])
kw = dict(ncol=1,
          fancybox=True, fontsize=10, frameon=False)
# leg1 = ax.legend(h[:], l[:], bbox_to_anchor=[0.5, 1.08], **kw)
ax.legend(handles=h, loc='lower right', **kw)

plt.axhline(1.675, c='tab:gray', lw=0.8)
plt.text(0.1, 1.375, r'$\delta_{\rm c}|_{G_{\rm eff}=1}$',
         fontsize=13, c='tab:grey')
ax.fill_between([-0.1, 1.1], 1.675-0.025, 1.675 +
                0.025, alpha=0.4, color='tab:gray')
plt.xlim(-0.05, 1.05)
"""



#plt.errorbar(x, y, yerr=yerr, xerr=xerr, fmt='.')
from multiprocessing import Pool
#############################################
#
# extract spectra and plot them
#
#############################################
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.integrate import odeint
from scipy import integrate
from scipy.special import lambertw
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.pylab as pl


plt.cla()
plt.figure()
plt.rcParams.update({"text.usetex": True})
fig = plt.figure(figsize=(4.25*2*.95, 2*2*1.1))


nn = 1
z_smf_arr = [0,1]

for z_smf in z_smf_arr:
    ax_Pk = plt.subplot(2,2,nn)

    
    ax_Pk.xaxis.set_minor_locator(AutoMinorLocator())
    ax_Pk.yaxis.set_minor_locator(AutoMinorLocator())


    plt.tick_params(axis='both', which='major', direction="in",
                    labelsize=14, length=5, top=True, right=True)
    plt.tick_params(axis='both', which='minor', direction="in",
                    labelsize=11, length=4, top=True, right=True)
    plt.tick_params(axis='both', which='major',
                    direction="in", labelsize=14, length=5)
    plt.tick_params(axis='both', which='minor',
                    direction="in", labelsize=11, length=4)



    obs = number_density(feature='GSMF', z_target=z_smf, h=h)
    j_data = 0
    k_func = 0
    colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                    '#ff7f00','#a65628','#f781bf','#999999']*4
    color_maps     = ['Reds', 'Blues', 'Greens'] *4
    markers        = ['o','s','v','^','<','>','p','*','D','.','8']*4
    linestyles     = ['-','--','-.',':']*4

    for ii in range(obs.n_target_observation):
        data       = obs.target_observation['Data'][ii]
        label      = obs.target_observation.index[ii]
        datatype   = obs.target_observation['DataType'][ii]
        color = 'tab:gray'
        marker     = 'o'
        linestyle  = linestyles[k_func]
        if datatype == 'data':
            if  ii == 0:
                ax_Pk.errorbar(10**data[:,0],  data[:,1],yerr = np.abs([data[:,1]-data[:,3],data[:,2]- data[:,1]]),\
                        label=r'$\rm pre-JWST$',capsize=0,ecolor=color,color='w',marker=marker,markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color)
            else:
                ax_Pk.errorbar(10**data[:,0],  data[:,1],yerr = np.abs([data[:,1]-data[:,3],data[:,2]- data[:,1]]),\
                        capsize=0,ecolor=color,color='w',marker=marker,markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color)
            
            j_data +=1
        
        color = 'k'
        ax_Pk.errorbar(1,1,yerr=1,label=r'$\rm JWST$',capsize=0,ecolor=color,color='w',marker='v',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color)

    if z_smf in [4,5,6,7,8]:
        path = '../observational_data/GSMF'
        Navarro = np.loadtxt(path + "/Navarro_z"+str(z_smf)+".dat")
        x = 10**Navarro[:,0]
        y = 1e-4*Navarro[:,1]
        yerr = 1e-4*Navarro[:,2]
        color = 'k'
        ax_Pk.errorbar(x,y,yerr=yerr,capsize=0,ecolor=color,color='w',marker='v',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color)
    
    pool_cpu = Pool(8)

    model = 'DES'
    model_H = 'LCDM'
    model_SFR = 'Puebla'

    pars2 = np.linspace(-1, 1, 10)
    pars1 = np.array([-0.75, 0, 1])
    n = len(pars2)
    cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#66c2a5"]) 
    cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#fc8d62"]) 
    cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#8da0cb"])

    colors = np.array([cmap1(np.linspace(0, 1, n)), cmap2(np.linspace(0, 1, n)), cmap3(np.linspace(0, 1, n))])
    f0 = 0.21
    
    for j, par1 in enumerate(pars1):
        SMF_library = SMF(1/(1+z_smf), model, model_H, model_SFR, pars1, pars2, 1e8, f0)
        Pk_arr = []
        for par2 in pars2:
            HMF_library = HMF(1/(1+z_smf), model, model_H, par1, par2, 1e8)
            Pk_arr.append(np.array(HMF_library.Pk(1/(1+z_smf), model, par1, par2))*h**3)
        k = kvec/h
        Masses = np.logspace(6,18,100)


        iterable = [(Masses, rhom, 1/(1+z_smf), model_H, model, model_SFR, par1, par2, k, Pk_arr[i], f0) for i,par2 in enumerate(pars2)]
        Masses_star, SMF_obs = zip(*pool_cpu.starmap(SMF_library.SMF_obs,tqdm(iterable, total=len(pars2))))
        for i in range(len(SMF_obs)):
            ax_Pk.plot(Masses_star[i], SMF_obs[i],c=colors[j][i])
    
    norm = plt.Normalize(pars2.min(), pars2.max())
    cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Grays, norm=norm), ax=ax_Pk)
    cbar.set_label(r'$T_2$', fontsize=16)

    #plt.errorbar(x.get('Duncan'),y.get('Duncan'),yerr=[yerr_down.get('Duncan'),yerr_up.get('Duncan')], c = 'tab:orange', capsize = 2, ls = 'None', marker = '.', label = r'$\rm Duncan+14$')
    #plt.errorbar(x.get('Song'),y.get('Song'),yerr=[yerr_down.get('Song'),yerr_up.get('Song')], c = 'tab:orange', capsize = 2, ls = 'None', marker = 's', label = r'$\rm Song+16$')
    #plines = plt.errorbar(x.get('Duncan'),y.get('Duncan'),yerr=[yerr_down.get('Duncan'),yerr_up.get('Duncan')],capsize=0,ecolor='tab:blue',color='w',marker='o',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='tab:blue')
    #plines = plt.errorbar(x.get('Song'),y.get('Song'),yerr=[yerr_down.get('Song'),yerr_up.get('Song')],capsize=0,ecolor='tab:orange',color='w',marker='s',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='tab:orange')


    #plines = plt.errorbar(x.get('Navarro'),y.get('Navarro'),yerr=[yerr_down.get('Navarro'),yerr_up.get('Navarro')],capsize=0,ecolor='k',color='w',marker=markers[j_data+1],markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='k', label = r'$\rm Navarro+2024$')

    # plt.scatter(1/a_vir-1, vir2, c = 'tab:orange')
    plt.xscale('log')
    plt.yscale('log')
    plt.tight_layout()
    plt.ylim(1e-8,1e0)
    plt.xlim(1e6,10**12.5)
    ax_Pk.set_xticklabels([])
    ax_Pk.set_ylabel(r'$\phi_{\star}\;[\rm Mpc^{-3}\;dex^{-1}]$', size = '16')

    plt.grid(".")
    
    ax_Pk.text(10**11,0.065,r'$z='+str(int(round(z_smf)))+r'$', size = '15')
    

    #hhh, llll = ax_Pk.get_legend_handles_labels()
    hhh = []

    line1 = Line2D([0], [0], label=r'$T_1=-1$', color='#398e73')
    line2 = Line2D([0], [0], label=r'$T_1=0$', color='#e64304')
    line3 = Line2D([0], [0], label=r'$T_1=1$', color='#48639e')
    hhh.extend([line1, line2, line3])
    kw = dict(ncol=1,
            fancybox=True, fontsize=10, frameon=False)
    # leg1 = ax.legend(h[:], l[:], bbox_to_anchor=[0.5, 1.08], **kw)
    ax_Pk.legend(handles=hhh, **kw,loc='lower left')



    nn += 1


z_smf_arr = [0,1]

for z_smf in z_smf_arr:
    ax_Pk = plt.subplot(2,2,nn)

    
    ax_Pk.xaxis.set_minor_locator(AutoMinorLocator())
    ax_Pk.yaxis.set_minor_locator(AutoMinorLocator())


    plt.tick_params(axis='both', which='major', direction="in",
                    labelsize=14, length=5, top=True, right=True)
    plt.tick_params(axis='both', which='minor', direction="in",
                    labelsize=11, length=4, top=True, right=True)
    plt.tick_params(axis='both', which='major',
                    direction="in", labelsize=14, length=5)
    plt.tick_params(axis='both', which='minor',
                    direction="in", labelsize=11, length=4)



    obs = number_density(feature='GSMF', z_target=z_smf, h=h)
    j_data = 0
    k_func = 0
    colors         = ['#e41a1c','#377eb8','#4daf4a','#984ea3',\
                    '#ff7f00','#a65628','#f781bf','#999999']*4
    color_maps     = ['Reds', 'Blues', 'Greens'] *4
    markers        = ['o','s','v','^','<','>','p','*','D','.','8']*4
    linestyles     = ['-','--','-.',':']*4

    for ii in range(obs.n_target_observation):
        data       = obs.target_observation['Data'][ii]
        label      = obs.target_observation.index[ii]
        datatype   = obs.target_observation['DataType'][ii]
        color = 'tab:gray'
        marker     = 'o'
        linestyle  = linestyles[k_func]
        if datatype == 'data':
            if  ii == 0:
                ax_Pk.errorbar(10**data[:,0],  data[:,1],yerr = np.abs([data[:,1]-data[:,3],data[:,2]- data[:,1]]),\
                        label=r'$\rm pre-JWST$',capsize=0,ecolor=color,color='w',marker=marker,markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color)
            else:
                ax_Pk.errorbar(10**data[:,0],  data[:,1],yerr = np.abs([data[:,1]-data[:,3],data[:,2]- data[:,1]]),\
                        capsize=0,ecolor=color,color='w',marker=marker,markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color)
            
            j_data +=1
        
        color = 'k'
        ax_Pk.errorbar(1,1,yerr=1,label=r'$\rm JWST$',capsize=0,ecolor=color,color='w',marker='v',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color)

    if z_smf in [4,5,6,7,8]:
        path = '../observational_data/GSMF'
        Navarro = np.loadtxt(path + "/Navarro_z"+str(z_smf)+".dat")
        x = 10**Navarro[:,0]
        y = 1e-4*Navarro[:,1]
        yerr = 1e-4*Navarro[:,2]
        color = 'k'
        ax_Pk.errorbar(x,y,yerr=yerr,capsize=0,ecolor=color,color='w',marker='v',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color)
    
    pool_cpu = Pool(8)

    model = 'DES'
    model_H = 'LCDM'
    model_SFR = 'double_power'

    pars2 = np.linspace(-1, 1, 10)
    pars1 = np.array([-0.75, 0, 1])
    n = len(pars2)
    cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#66c2a5"]) 
    cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#fc8d62"]) 
    cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#8da0cb"])

    colors = np.array([cmap1(np.linspace(0, 1, n)), cmap2(np.linspace(0, 1, n)), cmap3(np.linspace(0, 1, n))])
    f0 = 0.21
    
    for j, par1 in enumerate(pars1):
        SMF_library = SMF(1/(1+z_smf), model, model_H, model_SFR, pars1, pars2, 1e8, f0)
        Pk_arr = []
        for par2 in pars2:
            HMF_library = HMF(1/(1+z_smf), model, model_H, par1, par2, 1e8)
            Pk_arr.append(np.array(HMF_library.Pk(1/(1+z_smf), model, par1, par2))*h**3)
        k = kvec/h
        Masses = np.logspace(6,18,100)


        iterable = [(Masses, rhom, 1/(1+z_smf), model_H, model, model_SFR, par1, par2, k, Pk_arr[i], f0) for i,par2 in enumerate(pars2)]
        Masses_star, SMF_obs = zip(*pool_cpu.starmap(SMF_library.SMF_obs,tqdm(iterable, total=len(pars2))))
        for i in range(len(SMF_obs)):
            ax_Pk.plot(Masses_star[i], SMF_obs[i],c=colors[j][i])
    
    norm = plt.Normalize(pars2.min(), pars2.max())
    cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Grays, norm=norm), ax=ax_Pk)
    cbar.set_label(r'$T_2$', fontsize=16)

    #plt.errorbar(x.get('Duncan'),y.get('Duncan'),yerr=[yerr_down.get('Duncan'),yerr_up.get('Duncan')], c = 'tab:orange', capsize = 2, ls = 'None', marker = '.', label = r'$\rm Duncan+14$')
    #plt.errorbar(x.get('Song'),y.get('Song'),yerr=[yerr_down.get('Song'),yerr_up.get('Song')], c = 'tab:orange', capsize = 2, ls = 'None', marker = 's', label = r'$\rm Song+16$')
    #plines = plt.errorbar(x.get('Duncan'),y.get('Duncan'),yerr=[yerr_down.get('Duncan'),yerr_up.get('Duncan')],capsize=0,ecolor='tab:blue',color='w',marker='o',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='tab:blue')
    #plines = plt.errorbar(x.get('Song'),y.get('Song'),yerr=[yerr_down.get('Song'),yerr_up.get('Song')],capsize=0,ecolor='tab:orange',color='w',marker='s',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='tab:orange')


    #plines = plt.errorbar(x.get('Navarro'),y.get('Navarro'),yerr=[yerr_down.get('Navarro'),yerr_up.get('Navarro')],capsize=0,ecolor='k',color='w',marker=markers[j_data+1],markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='k', label = r'$\rm Navarro+2024$')

    # plt.scatter(1/a_vir-1, vir2, c = 'tab:orange')
    plt.xscale('log')
    plt.yscale('log')
    plt.tight_layout()
    plt.ylim(1e-8,1e0)
    plt.xlim(1e6,10**12.5)
    ax_Pk.set_xlabel(r'$M_\star\;[M_\odot]$', size = '16')
    ax_Pk.set_ylabel(r'$\phi_{\star}\;[\rm Mpc^{-3}\;dex^{-1}]$', size = '16')
    plt.grid(".")
    
    ax_Pk.text(10**11,0.065,r'$z='+str(int(round(z_smf)))+r'$', size = '15')
    



    #hhh, llll = ax_Pk.get_legend_handles_labels()
    hhh = []
    line1 = Line2D([0], [0], label=r'$T_1=-1$', color='#398e73')
    line2 = Line2D([0], [0], label=r'$T_1=0$', color='#e64304')
    line3 = Line2D([0], [0], label=r'$T_1=1$', color='#48639e')
    hhh.extend([line1, line2, line3])
    kw = dict(ncol=1,
            fancybox=True, fontsize=10, frameon=False)
    # leg1 = ax.legend(h[:], l[:], bbox_to_anchor=[0.5, 1.08], **kw)
    ax_Pk.legend(handles=hhh, **kw,loc='lower left')



    nn += 1


"""ac_arr = np.linspace(0.01, 1, 15)
par1 = 500
par2 = 0
deltac = delta_c(ac_arr, model, model_H, par1, par2)


di = deltac.binary_search_di(1, model, model_H, par1, par2,
                             0, len(delta_ini)-1, abs_err)
print(di)
delta_nl = deltac.collapse(di, model, model_H, par1, par2)
delta = delta_nl[:, 1]
a = delta_nl[:, 0]
cosmological_library = cosmological_functions(
    a, model, model_H, par1, par2)
H = cosmological_library.H_f(a, model_H, par1, par2)
dH = cosmological_library.dH_f(a, model_H, par1, par2)


G = 1/(8*np.pi)
Hdot = a*H*dH
beta = 1 + 2*H*par1/c*(1+Hdot/(3*H**2))
epsilon = 8/(9*beta**2)*(H0*par1/c)**2*Omegam0*a**(-3)
RRV = (epsilon*delta)**(-1/3)
mu = cosmological_library.mu(
    a, model, model_H, par1, par2, type='nonlinear', x=RRV)

plt.plot(a, epsilon*delta)


par1 = 3000
par2 = 0
deltac = delta_c(ac_arr, model, model_H, par1, par2)


di = deltac.binary_search_di(1, model, model_H, par1, par2,
                             0, len(delta_ini)-1, abs_err)
print(di)
delta_nl = deltac.collapse(di, model, model_H, par1, par2)
delta = delta_nl[:, 1]
a = delta_nl[:, 0]
cosmological_library = cosmological_functions(
    a, model, model_H, par1, par2)
H = cosmological_library.H_f(a, model_H, par1, par2)
dH = cosmological_library.dH_f(a, model_H, par1, par2)


G = 1/(8*np.pi)
Hdot = a*H*dH
beta = 1 + 2*H*par1/c*(1+Hdot/(3*H**2))
epsilon = 8/(9*beta**2)*(H0*par1/c)**2*Omegam0*a**(-3)
RRV = (epsilon*delta)**(-1/3)
mu = cosmological_library.mu(
    a, model, model_H, par1, par2, type='nonlinear', x=RRV)

plt.plot(a, epsilon*delta)
plt.xlim(0.3, 1)
plt.ylim(0.05, 500)
plt.yscale('log')

"""
"""model = 'nDGP'
model_H = 'nDGP'
par1 = 3000
par2 = 0.5
a_arr = np.linspace(0.1, 1, 10)
reion = reionization(a_arr, model, model_H, par1, par2)

R_arr = reion.radius_solve(model, model_H, par1, par2, a_arr)

plt.plot(a_arr, R_arr)
plt.yscale('log')"""
# plt.plot(a, mu)
"""

model = 'kmoufl'
model_H = 'kmoufl'
par2 = 0.5
ac_arr = np.array([0.1, 0.5, 0.75, 1])
pars1 = np.linspace(0.1, 0.3, 10)

n = len(ac_arr)
colors = pl.cm.Reds(np.linspace(0, 1, 5))


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
    plt.plot(pars1, Delta, c=colors[i+1])

plt.xscale('log')

plt.axhline(18*np.pi**2, ls=':', c='tab:blue')"""
"""
par1 = 3000
par2 = 0
ac_arr = np.linspace(0.01, 1, 15)

for ac in ac_arr:
    deltac = delta_c(ac, model, model_H, par1, par2)
    dc = deltac.delta_c_at_ac(ac, model, model_H, par1, par2)
    plt.scatter(ac, dc, c='tab:blue')

plt.axhline(1.688, c='tab:blue', ls=':')


par1 = 500
model_H = 'nDGP'
model_SFR = 'Puebla'

    deltac = delta_c(ac, model, model_H, par1, par2)
    dc = deltac.delta_c_at_ac(ac, model, model_H, par1, par2)
    plt.scatter(ac, dc, c='tab:orange')

plt.axhline(1.687, c='tab:orange', ls=':')
"""
"""
model = 'nDGP'
par2 = 0.5
ac_arr = np.array([0.1, 0.5, 0.75, 1])
pars1 = np.logspace(2.69897000434, 5, 10)

n = len(ac_arr)
colors = pl.cm.Reds(np.linspace(0, 1, 5))


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
    plt.plot(pars1, Delta, c=colors[i+1])

plt.xscale('log')

plt.axhline(18*np.pi**2, ls=':', c='tab:blue')"""
"""
par2 = 0.3
ac_arr = np.linspace(0.01, 1, 15)
par1 = 3000


deltac = delta_c(ac_arr, model, model_H, par1, par2)


di = deltac.binary_search_di(1, model, model_H, par1, par2,
                             0, len(delta_ini)-1, abs_err)
print(di)
delta_nl = deltac.collapse(di, model, model_H, par1, par2)
delta = delta_nl[:, 1]
a = delta_nl[:, 0]
cosmological_library = cosmological_functions(
    a, model, model_H, par1, par2)
H = cosmological_library.H_f(a, model_H, par1, par2)
dH = cosmological_library.dH_f(a, model_H, par1, par2)


G = 1/(8*np.pi)
Hdot = a*H*dH
beta = 1 + 2*H*par1/c*(1+Hdot/(3*H**2))
epsilon = 8/(9*beta**2)*(H0*par1/c)**2*Omegam0*a**(-3)
RRV = (epsilon*delta)**(-1/3)
mu = cosmological_library.mu(
    a, model, model_H, par1, par2, type='nonlinear', x=RRV)

plt.plot(a, mu, ls='-', c='tab:orange')
mu = cosmological_library.mu(
    a, model, model_H, par1, par2, type='linear', x=RRV)

plt.plot(a, mu, ls=':', c='tab:orange')


par2 = 0.3
ac_arr = np.linspace(0.01, 1, 15)
par1 = 500


deltac = delta_c(ac_arr, model, model_H, par1, par2)


di = deltac.binary_search_di(1, model, model_H, par1, par2,
                             0, len(delta_ini)-1, abs_err)
print(di)
delta_nl = deltac.collapse(di, model, model_H, par1, par2)
delta = delta_nl[:, 1]
a = delta_nl[:, 0]
cosmological_library = cosmological_functions(
    a, model, model_H, par1, par2)
H = cosmological_library.H_f(a, model_H, par1, par2)
dH = cosmological_library.dH_f(a, model_H, par1, par2)


G = 1/(8*np.pi)
Hdot = a*H*dH
beta = 1 + 2*H*par1/c*(1+Hdot/(3*H**2))
epsilon = 8/(9*beta**2)*(H0*par1/c)**2*Omegam0*a**(-3)
RRV = (epsilon*delta)**(-1/3)
mu = cosmological_library.mu(
    a, model, model_H, par1, par2, type='nonlinear', x=RRV)

plt.plot(a, mu, ls='-', c='tab:blue')
mu = cosmological_library.mu(
    a, model, model_H, par1, par2, type='linear', x=RRV)

plt.plot(a, mu, ls=':', c='tab:blue')


plt.yscale('log')
plt.xlim(0.2, 1)"""
# p
"""
plt.yscale('log')
plt.ylabel(r'$\delta_{\rm m}(a)$', size='16')
# plt.ylim(1.65, 1.7)
# plt.axhline(1.686, c='tab:gray')

# plt.axhline(1.672, c='tab:gray', ls=':')
# plt.xlim(10**(-3),1)
# plt.legend(loc='best')
plt.grid(".")
"""


plt.tight_layout()
plt.savefig('SMF_DES.pdf', bbox_inches='tight')
