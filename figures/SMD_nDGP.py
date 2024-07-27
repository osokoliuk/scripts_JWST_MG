

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

plt.cla()
plt.figure()
plt.rcParams.update({"text.usetex": True})
fig = plt.figure(figsize=(4.25*2*.95*0.9, 2*5*1.05*0.9))


nn = 1
z_smf_arr = [8,9,8,9]

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


    if z_smf == 8:
        f = np.genfromtxt("../observational_data/Desprez2024_MassiveGal.dat", names=True)
        sel = f['type'] == 0
        ax_Pk.errorbar(f['Mstar'][sel], f['Phi'][sel], 
                    marker = '+', linestyle='none', color='darkcyan', mec='darkcyan', ms=12, mew=2,
                    label=r"$\rm Desprez+24$")
        sel = f['type'] == 1
        ax_Pk.plot(f['Mstar'][sel], f['Phi'][sel], '-', color='darkcyan', mec='darkcyan', lw=3, alpha=0.6)

        ##
        f = np.genfromtxt("../observational_data/Labbe2023_MassiveGal.dat", names=True)
        sel = f['z'] == 8
        ax_Pk.errorbar(f['Mstar'][sel], f['Phi'][sel] + np.log10(factor_from_Salpeter_to_Chabrier), yerr=(f['Phi'][sel] - f['Philo'][sel], f['Phiup'][sel] - f['Phi'][sel]), 
                    marker = 's', linestyle='none', color='k', mec='k', capthick=2, capsize=6, ms=12)

        alpha = -1.82;  Phis = 10**(-4.69); Ms = 10**9.98
        mlim = 10**np.linspace(8, 12, 101)
        ax_Pk.plot(np.log10(mlim), np.log10(cumulative(mlim, Phis, Ms, alpha)) + np.log10(factor_from_Salpeter_to_Chabrier), 
                '-', lw=12, color='violet', alpha=0.4, zorder=-5, label=r"$\rm Stefanon+21$")

        alpha = -2.16;  Phis = 10**(-4.86); Ms = 10**10.0
        mlim = 10**np.linspace(8, 12, 101)
        ax_Pk.plot(np.log10(mlim), np.log10(cumulative(mlim, Phis, Ms, alpha)) + np.log10(0.63/0.67), 
                '-', lw=12, color='orange', alpha=0.4, zorder=-5)

        ##
        ax_Pk.errorbar(10.07681, 4.54678, yerr=([4.54678 - 3.32218], [5.07004 - 4.54678]), marker='o', 
                color='chocolate', mec='chocolate', ms=12, mew=2, lw=2, capsize=6, capthick=2, label=r"$\rm Wang+24$")

        ax_Pk.errorbar(10.85853, 5.32855, yerr=([5.32855 - 4.10212], [5.84263 - 5.32855]), marker='*', 
                color='chocolate', mec='chocolate', ms=16, mew=2, lw=2, capsize=6, capthick=2)

        ##
        f = np.genfromtxt("./observational_data/Akins2023_MassiveGal.dat", names=True)
        ax_Pk.errorbar(f['logMstar'], f['Phi'] + np.log10(0.63/0.67), yerr=(f['Phi'] - f['lo'], f['up'] - f['Phi']), 
                    xerr = (f['logMstar'] - f['left'], f['right'] - f['logMstar']),
                    marker = '^', linestyle='none', color='dimgray', mec='dimgray', capthick=2, capsize=6, ms=12, label=r"$\rm Akins+23$")
    else:
        data = np.loadtxt('/home/oleksii/codes/scripts_JWST_MG/observational_data/SMD/JWST_z9.txt')


    pool_cpu = Pool(8)

    model = 'nDGP'
    model_H = 'nDGP'
    if nn >= 3:
        model_SFR = 'double_power'
    else:
        model_SFR = 'Puebla'

    pars1 = np.logspace(2.5,4,10)
    par2 = 0
    f0 = 0.21
    n = len(pars1)
    cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#8da0cb"])

    colors = cmap3(np.linspace(0, 1, n))
    
    """    SMD_library = SMD(1/(1+z_smf), model, model_H, model_SFR, pars1, par2, 1e8, f0)
    Pk_arr = []
    for par1 in pars1:
        HMF_library = HMF(1/(1+z_smf), model, model_H, par1, par2, 1e8)
        Pk_arr.append(np.array(HMF_library.Pk(1/(1+z_smf), model, par1, par2))*h**3)
    k = kvec/h
    Masses = np.logspace(6,16,150)


    iterable = [(Masses, rhom, 1/(1+z_smf), model_H, model, model_SFR, par1, par2, k, Pk_arr[i], f0) for i,par1 in enumerate(pars1)]
    Masses_stars, SMDs = zip(*pool_cpu.starmap(SMD_library.SMD,tqdm(iterable, total=len(pars1))))

   #for i in range(len(SMF_obs)):
    
    for i in range(len(Masses_stars)):
        ax_Pk.plot(np.log10(Masses_stars[i]), np.log10(SMDs[i]), c = colors[i], lw=  1.5)
    """    #ax_Pk.fill_between(np.log10(Masses_star[2]), np.log10(SMF_obs[0]), np.log10(SMF_obs[2]), color='tab:gray', alpha=0.3)

    #plt.errorbar(x.get('Duncan'),y.get('Duncan'),yerr=[yerr_down.get('Duncan'),yerr_up.get('Duncan')], c = 'tab:orange', capsize = 2, ls = 'None', marker = '.', label = r'$\rm Duncan+14$')
    #plt.errorbar(x.get('Song'),y.get('Song'),yerr=[yerr_down.get('Song'),yerr_up.get('Song')], c = 'tab:orange', capsize = 2, ls = 'None', marker = 's', label = r'$\rm Song+16$')
    #plines = plt.errorbar(x.get('Duncan'),y.get('Duncan'),yerr=[yerr_down.get('Duncan'),yerr_up.get('Duncan')],capsize=0,ecolor='tab:blue',color='w',marker='o',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='tab:blue')
    #plines = plt.errorbar(x.get('Song'),y.get('Song'),yerr=[yerr_down.get('Song'),yerr_up.get('Song')],capsize=0,ecolor='tab:orange',color='w',marker='s',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='tab:orange')


    #plines = plt.errorbar(x.get('Navarro'),y.get('Navarro'),yerr=[yerr_down.get('Navarro'),yerr_up.get('Navarro')],capsize=0,ecolor='k',color='w',marker=markers[j_data+1],markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='k', label = r'$\rm Navarro+2024$')

    # plt.scatter(1/a_vir-1, vir2, c = 'tab:orange')
    plt.tight_layout()
    plt.xlim(5,11.5)
    plt.ylim(2.5,7.75)
    if nn != 3 and nn != 4:
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
            ax_Pk.set_ylabel(r'$\log_{10}\rho_\star\;[M_\odot\;\rm Mpc^{-3}]$', size='16')
    else:
        if nn % 2 == 0:
            ax_Pk.set_xlabel(r'$ \log_{10}M_\star\;[M_\odot]$', size = '16')
            ax_Pk.set_yticklabels([])
        else:
            ax_Pk.set_xlabel(r'$ \log_{10}M_\star\;[M_\odot]$', size = '16')
            ax_Pk.set_ylabel(r'$\log_{10}\rho_\star\;[M_\odot\;\rm Mpc^{-3}]$', size='16')

    plt.grid(".")
    
    ax_Pk.text(0.8,0.85,r'$z='+str(int(round(z_smf)))+r'$', size = '15', transform=ax_Pk.transAxes)
    
    legend1 = ax_Pk.legend(loc='lower left',fancybox=True, fontsize=10)
    legend1.get_frame().set_facecolor('none')
    legend1.get_frame().set_linewidth(0.0)
    ax_Pk.add_artist(legend1)
    
    nn += 1



mpl.rcParams['font.family'] = 'sans-serif'

#norm = colorss.Norm(pars1.min(), pars1.max())
norm = mpl.colors.Normalize(vmin=np.log10(pars1.min()), vmax=np.log10(pars1.max()))
ax_cbar = fig.add_axes([0.0805, 0.9835, 0.899, 0.01])
cbar_ax = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap3, norm=norm), cax=ax_cbar, orientation='horizontal', location = 'top', ticks=LinearLocator(numticks=8))
cbar_ax.set_label(r'$\log_{10}r_c$', fontsize=16)
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


plt.savefig('SMD_nDGP.pdf', bbox_inches='tight')
