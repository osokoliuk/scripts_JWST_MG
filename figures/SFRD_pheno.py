

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
from JWST_MG.SMD import SMD
from JWST_MG.SMF import SMF
from JWST_MG.UVLF import UVLF
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
fig = plt.figure(figsize=(4.25*1*.95*0.9, 1*2*0.95*0.9))
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
cbar.set_label(r'$E_{11}$', fontsize=16)

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
cbar = plt.colorbar(mpfigures/SMF_E11.py
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
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import LinearLocator
import matplotlib.pylab as pl
mpl.rcParams['axes.linewidth'] = 1.5

plt.clf()
plt.figure()
plt.rcParams.update({"text.usetex": True})
fig = plt.figure(figsize=(4.25*2*.95*0.9, 2*5*1.05*0.9))


nn = 1
z_smf_arr = [0,1,2,3,0,1,2,3]
pool_cpu = Pool(8)


def FindValueIndex(seq, val):
    r = np.where(np.diff(np.sign(seq - val)) != 0)
    idx = r + (val - seq[r]) / (seq[r + np.ones_like(r)] - seq[r])
    idx = np.append(idx, np.where(seq == val))
    idx = np.sort(idx)
    return int(np.round(idx)[0])



for z_smf in z_smf_arr:
    ax_Pk = plt.subplot(4,2,nn)

    
    ax_Pk.xaxis.set_minor_locator(AutoMinorLocator())
    ax_Pk.yaxis.set_minor_locator(AutoMinorLocator())


    plt.tick_params(axis='both', which='major', direction="in",
                    labelsize=14, length=5, top=True, right=True, width = 1.5)
    plt.tick_params(axis='both', which='minor', direction="in",
                    labelsize=11, length=4, top=True, right=True, width = 1.1)
    plt.tick_params(axis='both', which='major',
                    direction="in", labelsize=14, length=5, width = 1.5)
    plt.tick_params(axis='both', which='minor',
                    direction="in", labelsize=11, length=4, width = 1.1)


    #data = np.loadtxt('/home/oleksii/codes/scripts_JWST_MG/observational_data/SFRD/sfrd_MD14.txt')
    #x = data[:,0]
    #xerr = data[:,1]
    #y = data[:,2]
    #yerr = data[:,3]
    
    data = np.loadtxt('/home/oleksii/codes/scripts_JWST_MG/observational_data/SFRD/csfrs_new.dat')
    x = data[:,0]
    y = data[:,1]
    yerr = [data[:,2], data[:,3]]
    xerr = 0
    marker = '.'
    color = 'tab:gray'
    plt.errorbar(x,y,yerr=yerr,xerr=xerr, capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 3, label = r'$\rm pre-JWST$')


    data = np.loadtxt('/home/oleksii/codes/scripts_JWST_MG/observational_data/SFRD/sfrd_jwst.txt')
    x = data[:,0]
    y = data[:,1]
    yerr = [data[:,3], data[:,2]]

    color = 'k'
    marker='^'

    plt.errorbar(x,y,yerr=yerr, capsize=2,ecolor=color,color='w',marker='^',markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 3, label = r'$\rm JWST$')




    if nn % 2 == 0:
        model_SFR = 'double_power'
    else:
        model_SFR = 'Puebla'

    if nn >= 1 and nn <= 2:
        model = 'E11'
        model_H = 'LCDM'
        pars2 = np.linspace(-1,1,5)
        pars1 = np.array([-1, 0, 1])
        n = len(pars2)
        f0 = 0.21
    elif nn >= 3 and nn <= 4:
        model = 'gmu'
        model_H = 'LCDM'
        pars2 = np.linspace(-1,1,5)
        pars1 = np.array([-1, 0, 1])
        n = len(pars2)
        f0 = 0.21
    elif nn >= 5 and nn <= 6:
        model = 'DES'
        model_H = 'LCDM'
        pars2 = np.linspace(-0.75,1,5)
        pars1 = np.array([-1, 0, 1])
        n = len(pars2)
        f0 = 0.21
    elif nn >= 7 and nn <= 8:
        model = 'wCDM'
        model_H = 'wCDM'
        pars2 = np.linspace(0.4,0.6,5)
        pars1 = np.array([-1.5, -1, -0.5])
        n = len(pars2)
        f0 = 0.21
    else:
        raise Exception("Beyond subplot limit")


    cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#66c2a5"]) 
    cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#fc8d62"]) 
    cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#8da0cb"])

    cmap = np.array([cmap1, cmap2, cmap3])
    colors = [None]*3
    colors[0] = cmap[0]((np.linspace(0, 1, n)))
    colors[1] = cmap[1]((np.linspace(0, 1, n)))
    colors[2] = cmap[2]((np.linspace(0, 1, n)))

    for ll, par1 in enumerate(pars1):
        z_int = np.linspace(0,18,19) #np.linspace(12,5,35)
        a_arr = 1/(1+z_int)

        UVLF_library = UVLF(1/(1+z_int), model, model_H, model_SFR, pars1, pars2, 1e8, f0)
        Pk_arr = np.empty(shape=(len(pars2), len(z_int)), dtype='object')
        Masses_arr = np.empty(shape=(len(pars2), len(z_int)), dtype='object')
        for j, par2 in enumerate(pars2):
            for i, z_i in enumerate(z_int):
                HMF_library = HMF(1/(1+z_i), model, model_H, par1, par2, 1e8)
                Pk_arr[j][i] = np.array(HMF_library.Pk(1/(1+z_i), model, par1, par2))*h**3
                Masses = np.logspace(8,18,350)                
                SMF_library = SMF(1/(1+z_i), model, model_H, model_SFR, par1, par2, Masses, f0)
                Mstar = Omegab0/Omegam0*Masses*SMF_library.epsilon(Masses, model_SFR, 1/(1+z_i), f0)
                idx = FindValueIndex(Mstar,1e8)
                Masses_arr[j][i] = Masses[idx:-1]
        k = kvec/h


        iterable = [(a_arr, rhom, model, model_H, model_SFR, par1, par2, Masses_arr[i], k, Pk_arr[i], f0) for i,par2 in enumerate(pars2)]
        SFRD_obs = pool_cpu.starmap( UVLF_library.SFRD,tqdm(iterable, total=len(pars2)))
        data = np.savez('./data_folder/SFRD_'+str(model)+'_'+str(model_SFR) +'_'+str(par1)+'.npz', name1 = z_int, name2 = SFRD_obs)


        #print(SMF)
        #Masses_star = SMF[0]
        #SMF_obs = SMF[1]
        for i in range(len(SFRD_obs)):
            plt.plot(z_int, np.log10(SFRD_obs[i]), c = colors[ll][i], lw = 1)
    

    if nn >= 1 and nn <= 2:
        line3 = ax_Pk.plot([0], [0], label=r'$E_{11}=-1$', color='#66c2a5')   
        line2 = ax_Pk.plot([0], [0], label=r'$E_{11}=0$', color='#fc8d62')
        line1 = ax_Pk.plot([0], [0], label=r'$E_{11}=1$', color='#8da0cb')
    elif nn >= 3 and nn <= 4:
        line3 = ax_Pk.plot([0], [0], label=r'$g_{\mu}=-1$', color='#66c2a5')   
        line2 = ax_Pk.plot([0], [0], label=r'$g_{\mu}=0$', color='#fc8d62')
        line1 = ax_Pk.plot([0], [0], label=r'$g_{\mu}=1$', color='#8da0cb')
    elif nn >= 5 and nn <= 6:
        line3 = ax_Pk.plot([0], [0], label=r'$T_{1}=-1$', color='#66c2a5')   
        line2 = ax_Pk.plot([0], [0], label=r'$T_{1}=0$', color='#fc8d62')
        line1 = ax_Pk.plot([0], [0], label=r'$T_{1}=1$', color='#8da0cb')
    elif nn >= 7 and nn <= 8:
        line3 = ax_Pk.plot([0], [0], label=r'$w_{\Lambda}=-1.5$', color='#66c2a5')   
        line2 = ax_Pk.plot([0], [0], label=r'$w_{\Lambda}=-1$', color='#fc8d62')
        line1 = ax_Pk.plot([0], [0], label=r'$w_{\Lambda}=-0.5$', color='#8da0cb')

    legend1 = ax_Pk.legend(loc='upper right',fancybox=True, fontsize=10)
    legend1.get_frame().set_facecolor('none')
    legend1.get_frame().set_linewidth(0.0)
    ax_Pk.add_artist(legend1)
    

    #plt.errorbar(x.get('Duncan'),y.get('Duncan'),yerr=[yerr_down.get('Duncan'),yerr_up.get('Duncan')], c = 'tab:orange', capsize = 2, ls = 'None', marker = '.', label = r'$\rm Duncan+14$')
    #plt.errorbar(x.get('Song'),y.get('Song'),yerr=[yerr_down.get('Song'),yerr_up.get('Song')], c = 'tab:orange', capsize = 2, ls = 'None', marker = 's', label = r'$\rm Song+16$')
    #plines = plt.errorbar(x.get('Duncan'),y.get('Duncan'),yerr=[yerr_down.get('Duncan'),yerr_up.get('Duncan')],capsize=0,ecolor='tab:blue',color='w',marker='o',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='tab:blue')
    #plines = plt.errorbar(x.get('Song'),y.get('Song'),yerr=[yerr_down.get('Song'),yerr_up.get('Song')],capsize=0,ecolor='tab:orange',color='w',marker='s',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='tab:orange')


    #plines = plt.errorbar(x.get('Navarro'),y.get('Navarro'),yerr=[yerr_down.get('Navarro'),yerr_up.get('Navarro')],capsize=0,ecolor='k',color='w',marker=markers[j_data+1],markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='k', label = r'$\rm Navarro+2024$')

    # plt.scatter(1/a_vir-1, vir2, c = 'tab:orange')
    plt.xlim(-1,19)
    plt.ylim(-4.5,0)
    if nn != 7 and nn != 8:
        if nn % 2 == 0:
            ax_Pk.set_xticklabels([])
            ax_Pk.set_yticklabels([])
            if nn == 2:
                #cbar.ax.tick_params(size=8, width=2, direction='in')
                """ax_divider = make_axes_locatable(ax_Pk)
                norm = colorss.LogNorm(pars1.min(), pars1.max())
                cax = ax_divider.append_axes("top", size="7%", pad="2%")
                cb = fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap3, norm=norm), cax=cax, orientation="horizontal", location = 'top')
                cb.set_label(r'$r_c$', fontsize=16)
                fig.colorbar(mappable, ax=axs)"""
        else:
            ax_Pk.set_xticklabels([])
            ax_Pk.set_ylabel(r'$\log_{10}\rho_{\rm SFR}\;[M_\odot\;\rm Mpc^{-3}]$', size='16')
    else:
        if nn % 2 == 0:
            ax_Pk.set_xlabel(r'$ \mathrm{Redshift}\;z$', size = '16')
            ax_Pk.set_yticklabels([])
        else:
            ax_Pk.set_xlabel(r'$ \mathrm{Redshift}\;z$', size = '16')
            ax_Pk.set_ylabel(r'$\log_{10}\rho_{\rm SFR}\;[M_\odot\;\rm Mpc^{-3}]$', size='16')

    plt.grid(".")
    
    
    
    nn += 1

mpl.rcParams['font.family'] = 'sans-serif'

#norm = colorss.Norm(pars1.min(), pars1.max())
norm = mpl.colors.Normalize(vmin=-1, vmax = 1)
ax_cbar = fig.add_axes([0.101, 0.9775, 0.8785, 0.01])
cbar_ax = plt.colorbar(mpl.cm.ScalarMappable(cmap=mpl.cm.Greys, norm=norm), cax=ax_cbar, orientation='horizontal', location = 'top', ticks=LinearLocator(numticks=8))
cbar_ax.set_label(r'$E_{22}/g_\gamma/T_2/\gamma$', fontsize=16)
cbar_ax.ax.tick_params(width=1.5, length=5, which = 'major')
cbar_ax.ax.tick_params(width=1.1, length=4, which = 'minor')
cbar_ax.ax.xaxis.set_minor_locator(AutoMinorLocator())
cbar_ax.ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(r'$'+'%.2g'+r'$'))




plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=-0.0)



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


plt.savefig('SFRD_pheno.pdf', bbox_inches='tight')
