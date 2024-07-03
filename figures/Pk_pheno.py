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
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator,LinearLocator)
import matplotlib.pylab as pl
import matplotlib.ticker as mticker
# import necessary modules
# uncomment to get plots displayed in notebook
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
import math
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 1.5

cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#8da0cb"])
marker = '.'
plt.figure()
plt.rcParams.update({"text.usetex":True})
fig = plt.figure(figsize=(4.25*2*.95*0.9, 2*5*1.05*0.9))
ax_Pk = plt.subplot(421)
n = 10
cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#66c2a5"]) 
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#fc8d62"]) 
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#8da0cb"])

cmap = np.array([cmap1, cmap2, cmap3])
colors = [None]*3
colors[0] = cmap[0]((np.linspace(0, 1, n)))
colors[1] = cmap[1]((np.linspace(0, 1, n)))
colors[2] = cmap[2]((np.linspace(0, 1, n)))


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


############################################
#
# Varying parameter (others fixed to default)
#
var_name = 'N_ur'
var_array = [3.044]
var_num = len(var_array)
var_legend = r'$N_\mathrm{eff}$'
var_figname = 'neff'
#
# Constraints to be matched
#
# As explained in the "Neutrino cosmology" book, CUP, Lesgourgues et al., section 5.3, the goal is to vary
# - omega_cdm by a factor alpha = (1 + coeff*Neff)/(1 + coeff*3.046)
# - h by a factor sqrt*(alpha)
# in order to keep a fixed z_equality(R/M) and z_equality(M/Lambda)
#
omega_b = 0.0223828
omega_cdm_standard = 0.1201075
h_standard = 0.67810
#
# coefficient such that omega_r = omega_gamma (1 + coeff*Neff),
# i.e. such that omega_ur = omega_gamma * coeff * Neff:
# coeff = omega_ur/omega_gamma/Neff_standard 
# We could extract omega_ur and omega_gamma on-the-fly within th script, 
# but for simplicity we did a preliminary interactive run with background_verbose=2
# and we copied the values given in the budget output.
#
coeff = 1.70961e-05/2.47298e-05/3.044
print ("coeff=",coeff)
#
#############################################
#
# Fixed settings
#
############################################
#
# Varying parameter (others fixed to default)
#
var_name = 'N_ur'
var_array = [3.044]
var_num = len(var_array)
var_legend = r'$N_\mathrm{eff}$'
var_figname = 'neff'
#
# Constraints to be matched
#
# As explained in the "Neutrino cosmology" book, CUP, Lesgourgues et al., section 5.3, the goal is to vary
# - omega_cdm by a factor alpha = (1 + coeff*Neff)/(1 + coeff*3.046)
# - h by a factor sqrt*(alpha)
# in order to keep a fixed z_equality(R/M) and z_equality(M/Lambda)
#
omega_b = 0.0223828
omega_cdm_standard = 0.1201075
h_standard = 0.67810
#
# coefficient such that omega_r = omega_gamma (1 + coeff*Neff),
# i.e. such that omega_ur = omega_gamma * coeff * Neff:
# coeff = omega_ur/omega_gamma/Neff_standard 
# We could extract omega_ur and omega_gamma on-the-fly within th script, 
# but for simplicity we did a preliminary interactive run with background_verbose=2
# and we copied the values given in the budget output.
#
coeff = 1.70961e-05/2.47298e-05/3.044
print ("coeff=",coeff)
#
#############################################
#
# Fixed settings
#

common_settings = {'A_s':2.101e-9,
          'n_s':0.9665,
          'tau_reio':0.0561,
          'omega_b':0.02242,
          'omega_cdm':0.11933,
          'h':0.6766,
          'YHe':0.2425,
          'T_cmb':2.7255,
          'gauge':'newtonian', #FOR MGCLASS TO WORK, GAUGE NEEDS TO BE NEWTONIAN
          'k_pivot': 0.05,
          'mg_z_init': 10.000,
          'l_logstep': 1.025,
          'l_linstep':15,
          'P_k_max_1/Mpc':113.0,
          'l_switch_limber':9,
          'perturb_sampling_stepsize': 0.05,
          'output':'tCl,pCl,lCl,mPk',
          'l_max_scalars': 3000,
          'lensing': 'yes',
          'mg_ansatz':'plk_late'}

#
##############################################
#
# loop over varying parameter values
#
beta_arr = [-1,0,1]
K0_arr = np.linspace(-1,1,10)


M = {}
#
for j, beta in enumerate(beta_arr):
    for i, K0 in enumerate(K0_arr):
        common_settings['mg_E11'] = beta
        common_settings['mg_E22'] = K0
        #
        # rescale omega_cdm and h
        #
        #
        # call CLASS
        #
        M = Class()
        M.set(common_settings)

        M.compute()
        
        


        kvec = np.logspace(-4,np.log10(113),1000) # array of kvec in h/Mpc
        twopi = 2.*math.pi
        #
        # Create figures
        #

        #
        # loop over varying parameter values
        #
        ll = {}
        clM = {}
        clTT = {}
        pkM = {}
        legarray = []
        Tcmb = 2.72e6

        #
        #
        # deal with colors and legends
        #
        h = 0.6766
        #
        # get Cls
        #
        # store P(k) for common k values
        #
        pkM = []
        # The function .pk(k,z) wants k in 1/Mpc so we must convert kvec for each case with the right h 
        khvec = kvec*h # This is k in 1/Mpc
        for kh in khvec:
            pkM.append(M.pk(kh,0.)*h**3) 
        #    
        # plot P(k)
        #
        f = mticker.ScalarFormatter(useMathText=True)
        f.set_powerlimits((-6,6))
        ax_Pk.plot(kvec,pkM,
                        color=colors[j][i],#alpha=var_alpha,
                        linestyle='-')

    #
    # plot C_l^TT
    #
    
x = np.loadtxt('wmap.txt')[:,0]
y = np.loadtxt('wmap.txt')[:,1]
z = np.loadtxt('wmap.txt')[:,2]
plines_full = []
color = 'k'
plines =ax_Pk.errorbar(x,y,yerr=(z-y)*2,capsize=2,ecolor=color,color='w',marker='^',markersize=4,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 3)

plines_full.append(plines)

x = np.loadtxt('lya.txt')[:,0]
y = np.loadtxt('lya.txt')[:,1]
z = np.loadtxt('lya.txt')[:,2]

color = 'tab:gray'
plines = ax_Pk.errorbar(x,y,yerr=(z-y)*2,capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 3)
plines_full.append(plines)

legend1 = ax_Pk.legend(plines_full, [r"$\rm  WMAP/ACT$",  r"$\rm \mathrm{Ly-}\alpha\rm \,forest$"], loc='lower left',fancybox=True, fontsize=11)
legend1.get_frame().set_facecolor('none')
legend1.get_frame().set_linewidth(0.0)
ax_Pk.add_artist(legend1)

ax_Pk.set_ylabel(r'$P_{\rm m}(k)\;[(h^{-1}\mathrm{Mpc})^3]$', size = '16')


ax_Pk.set_xscale('log')
ax_Pk.set_yscale('log')
#plt.ylim(10**1.6,10**4.75)
#plt.xlim(10**(-3),1)
#plt.legend(loc='best')
ax_Pk.grid(".")
"""
h, l = ax_Pk.get_legend_handles_labels()
kw = dict(ncol=2, loc="lower center",fancybox=True, fontsize=11,frameon=False)    
leg2 = ax_Pk.legend(h[:],l[:], bbox_to_anchor=[0.5,1.08],**kw)
"""
ax_Pk.set_xscale('log')
ax_Pk.set_yscale('log')

ax_Pk.set_xlim(10**(-4),4*10**0)
ax_Pk.set_ylim(10**(0),5*10**5)

#plt.legend(loc='best')
plt.grid(".")
line3 = ax_Pk.plot([0], [0], label=r'$E_{11}=-1$', color='#66c2a5')   
line2 = ax_Pk.plot([0], [0], label=r'$E_{11}=0$', color='#fc8d62')
line1 = ax_Pk.plot([0], [0], label=r'$E_{11}=1$', color='#8da0cb')

legend1 = ax_Pk.legend(loc='upper right',fancybox=True, fontsize=10)
legend1.get_frame().set_facecolor('none')
legend1.get_frame().set_linewidth(0.0)
ax_Pk.add_artist(legend1)


pl.rcParams['font.family'] = 'sans-serif'
p2 = ax_Pk.get_position().get_points().flatten()
#norm = colorss.Norm(pars1.min(), pars1.max())
K0_arr = np.array(K0_arr)
norm = mpl.colors.Normalize(vmin=K0_arr.min(), vmax=K0_arr.max())
ax_cbar = fig.add_axes([p2[0]-0.0218, 0.9845, p2[2]-p2[0]+0.085, 0.01])
cbar_ax = plt.colorbar(mpl.cm.ScalarMappable(cmap=mpl.cm.Greys, norm=norm), cax=ax_cbar, orientation='horizontal', location = 'top', ticks=LinearLocator(numticks=8))
cbar_ax.set_label(r'$E_{22}$', fontsize=16)
cbar_ax.ax.tick_params(width=1.5, length=5, which = 'major')
cbar_ax.ax.tick_params(width=1.1, length=4, which = 'minor')
cbar_ax.ax.xaxis.set_minor_locator(AutoMinorLocator())
cbar_ax.ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(r'$'+'%.2g'+r'$'))
labels = cbar_ax.ax.get_xticklabels()
labels[0] = ""
cbar_ax.ax.set_xticklabels(labels)
ax_Pk.xaxis.set_ticklabels([])

ax_Pk = plt.subplot(422)
n = 10
cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#66c2a5"]) 
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#fc8d62"]) 
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#8da0cb"])

cmap = np.array([cmap1, cmap2, cmap3])
colors = [None]*3
colors[0] = cmap[0]((np.linspace(0, 1, n)))
colors[1] = cmap[1]((np.linspace(0, 1, n)))
colors[2] = cmap[2]((np.linspace(0, 1, n)))


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


############################################
#
# Varying parameter (others fixed to default)
#
var_name = 'N_ur'
var_array = [3.044]
var_num = len(var_array)
var_legend = r'$N_\mathrm{eff}$'
var_figname = 'neff'
#
# Constraints to be matched
#
# As explained in the "Neutrino cosmology" book, CUP, Lesgourgues et al., section 5.3, the goal is to vary
# - omega_cdm by a factor alpha = (1 + coeff*Neff)/(1 + coeff*3.046)
# - h by a factor sqrt*(alpha)
# in order to keep a fixed z_equality(R/M) and z_equality(M/Lambda)
#
omega_b = 0.0223828
omega_cdm_standard = 0.1201075
h_standard = 0.67810
#
# coefficient such that omega_r = omega_gamma (1 + coeff*Neff),
# i.e. such that omega_ur = omega_gamma * coeff * Neff:
# coeff = omega_ur/omega_gamma/Neff_standard 
# We could extract omega_ur and omega_gamma on-the-fly within th script, 
# but for simplicity we did a preliminary interactive run with background_verbose=2
# and we copied the values given in the budget output.
#
coeff = 1.70961e-05/2.47298e-05/3.044
print ("coeff=",coeff)
#
#############################################
#
# Fixed settings
#
############################################
#
# Varying parameter (others fixed to default)
#
var_name = 'N_ur'
var_array = [3.044]
var_num = len(var_array)
var_legend = r'$N_\mathrm{eff}$'
var_figname = 'neff'
#
# Constraints to be matched
#
# As explained in the "Neutrino cosmology" book, CUP, Lesgourgues et al., section 5.3, the goal is to vary
# - omega_cdm by a factor alpha = (1 + coeff*Neff)/(1 + coeff*3.046)
# - h by a factor sqrt*(alpha)
# in order to keep a fixed z_equality(R/M) and z_equality(M/Lambda)
#
omega_b = 0.0223828
omega_cdm_standard = 0.1201075
h_standard = 0.67810
#
# coefficient such that omega_r = omega_gamma (1 + coeff*Neff),
# i.e. such that omega_ur = omega_gamma * coeff * Neff:
# coeff = omega_ur/omega_gamma/Neff_standard 
# We could extract omega_ur and omega_gamma on-the-fly within th script, 
# but for simplicity we did a preliminary interactive run with background_verbose=2
# and we copied the values given in the budget output.
#
coeff = 1.70961e-05/2.47298e-05/3.044
print ("coeff=",coeff)
#
#############################################
#
# Fixed settings
#

common_settings = {'A_s':2.101e-9,
          'n_s':0.9665,
          'tau_reio':0.0561,
          'omega_b':0.02242,
          'omega_cdm':0.11933,
          'h':0.6766,
          'YHe':0.2425,
          'T_cmb':2.7255,
          'gauge':'newtonian', #FOR MGCLASS TO WORK, GAUGE NEEDS TO BE NEWTONIAN
          'k_pivot': 0.05,
          'mg_z_init': 10.000,
          'l_logstep': 1.025,
          'l_linstep':15,
          'P_k_max_1/Mpc':113.0,
          'l_switch_limber':9,
          'perturb_sampling_stepsize': 0.05,
          'output':'tCl,pCl,lCl,mPk',
          'l_max_scalars': 3000,
          'lensing': 'yes',
          'mg_ansatz':'z_flex_late'}

#
##############################################
#
# loop over varying parameter values
#
beta_arr = [-1,0,1]
K0_arr = np.linspace(-1,1,10)


labels = [r'$g_{\mu}=0$   $g_{\gamma}=0$',
          r'$g_{\mu}=0.1$ $g_{\gamma}=0.1$',
          r'$g_{\mu}=0.1$ $g_{\gamma}=0.2$',
          r'$g_{\mu}=0.1$ $g_{\gamma}=0.5$',
          r'$g_{\mu}=0.2$ $g_{\gamma}=0.1$',
          r'$g_{\mu}=0.2$ $g_{\gamma}=0.2$',
          r'$g_{\mu}=0.2$ $g_{\gamma}=0.5$',
          r'$g_{\mu}=0.5$ $g_{\gamma}=0.1$',
          r'$g_{\mu}=0.5$ $g_{\gamma}=0.2$',
          r'$g_{\mu}=0.5$ $g_{\gamma}=0.5$']

M = {}
#
for j, beta in enumerate(beta_arr):
    for i, K0 in enumerate(K0_arr):
        common_settings['mg_muz'] = beta
        common_settings['mg_gamz'] = K0
        common_settings['mg_zzn'] = 1

        #
        # rescale omega_cdm and h
        #
        #
        # call CLASS
        #
        M = Class()
        M.set(common_settings)

        M.compute()
        
        


        kvec = np.logspace(-4,np.log10(113),1000) # array of kvec in h/Mpc
        twopi = 2.*math.pi
        #
        # Create figures
        #

        #
        # loop over varying parameter values
        #
        ll = {}
        clM = {}
        clTT = {}
        pkM = {}
        legarray = []
        Tcmb = 2.72e6

        #
        #
        # deal with colors and legends
        #
        h = 0.6766
        #
        # get Cls
        #
        # store P(k) for common k values
        #
        pkM = []
        # The function .pk(k,z) wants k in 1/Mpc so we must convert kvec for each case with the right h 
        khvec = kvec*h # This is k in 1/Mpc
        for kh in khvec:
            pkM.append(M.pk(kh,0.)*h**3) 
        #    
        # plot P(k)
        #
        f = mticker.ScalarFormatter(useMathText=True)
        f.set_powerlimits((-6,6))
        ax_Pk.plot(kvec,pkM,
                        color=colors[j][i],#alpha=var_alpha,
                        linestyle='-')

    #
    # plot C_l^TT
    #
    
x = np.loadtxt('wmap.txt')[:,0]
y = np.loadtxt('wmap.txt')[:,1]
z = np.loadtxt('wmap.txt')[:,2]
plines_full = []
color = 'k'
plines =ax_Pk.errorbar(x,y,yerr=(z-y)*2,capsize=2,ecolor=color,color='w',marker='^',markersize=4,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 3)

plines_full.append(plines)

x = np.loadtxt('lya.txt')[:,0]
y = np.loadtxt('lya.txt')[:,1]
z = np.loadtxt('lya.txt')[:,2]

color = 'tab:gray'
plines = ax_Pk.errorbar(x,y,yerr=(z-y)*2,capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 3)
plines_full.append(plines)

legend1 = ax_Pk.legend(plines_full, [r"$\rm  WMAP/ACT$",  r"$\rm \mathrm{Ly-}\alpha\rm \,forest$"], loc='lower left',fancybox=True, fontsize=11)
legend1.get_frame().set_facecolor('none')
legend1.get_frame().set_linewidth(0.0)
ax_Pk.add_artist(legend1)

ax_Pk.set_xscale('log')
ax_Pk.set_yscale('log')
#plt.ylim(10**1.6,10**4.75)
#plt.xlim(10**(-3),1)
#plt.legend(loc='best')
ax_Pk.grid(".")
"""
h, l = ax_Pk.get_legend_handles_labels()
kw = dict(ncol=2, loc="lower center",fancybox=True, fontsize=11,frameon=False)    
leg2 = ax_Pk.legend(h[:],l[:], bbox_to_anchor=[0.5,1.08],**kw)
"""
ax_Pk.set_xscale('log')
ax_Pk.set_yscale('log')

ax_Pk.set_xlim(10**(-4),4*10**0)
ax_Pk.set_ylim(10**(0),5*10**5)
ax_Pk.yaxis.set_ticklabels([])

#plt.legend(loc='best')
plt.grid(".")
line3 = ax_Pk.plot([0], [0], label=r'$g_{\mu}=-1$', color='#66c2a5')   
line2 = ax_Pk.plot([0], [0], label=r'$g_{\mu}=0$', color='#fc8d62')
line1 = ax_Pk.plot([0], [0], label=r'$g_{\mu}=1$', color='#8da0cb')

legend1 = ax_Pk.legend(loc='upper right',fancybox=True, fontsize=10)
legend1.get_frame().set_facecolor('none')
legend1.get_frame().set_linewidth(0.0)
ax_Pk.add_artist(legend1)


pl.rcParams['font.family'] = 'sans-serif'
p2 = ax_Pk.get_position().get_points().flatten()
#norm = colorss.Norm(pars1.min(), pars1.max())
K0_arr = np.array(K0_arr)
norm = mpl.colors.Normalize(vmin=K0_arr.min(), vmax=K0_arr.max())
x = 0.0167
ax_cbar = fig.add_axes([p2[0]-0.0218+x, 0.9845, p2[2]-p2[0]+0.085, 0.01])
cbar_ax = plt.colorbar(mpl.cm.ScalarMappable(cmap=mpl.cm.Greys, norm=norm), cax=ax_cbar, orientation='horizontal', location = 'top', ticks=LinearLocator(numticks=8))
cbar_ax.set_label(r'$K_0$', fontsize=16)
cbar_ax.ax.tick_params(width=1.5, length=5, which = 'major')
cbar_ax.ax.tick_params(width=1.1, length=4, which = 'minor')
cbar_ax.ax.xaxis.set_minor_locator(AutoMinorLocator())
cbar_ax.ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(r'$'+'%.2g'+r'$'))
labels = cbar_ax.ax.get_xticklabels()
labels[0] = ""
cbar_ax.ax.set_xticklabels(labels)
ax_Pk.xaxis.set_ticklabels([])


ax_Pk = plt.subplot(423)
n = 10
cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#66c2a5"]) 
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#fc8d62"]) 
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#8da0cb"])

cmap = np.array([cmap1, cmap2, cmap3])
colors = [None]*3
colors[0] = cmap[0]((np.linspace(0, 1, n)))
colors[1] = cmap[1]((np.linspace(0, 1, n)))
colors[2] = cmap[2]((np.linspace(0, 1, n)))


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


############################################
#
# Varying parameter (others fixed to default)
#
var_name = 'N_ur'
var_array = [3.044]
var_num = len(var_array)
var_legend = r'$N_\mathrm{eff}$'
var_figname = 'neff'
#
# Constraints to be matched
#
# As explained in the "Neutrino cosmology" book, CUP, Lesgourgues et al., section 5.3, the goal is to vary
# - omega_cdm by a factor alpha = (1 + coeff*Neff)/(1 + coeff*3.046)
# - h by a factor sqrt*(alpha)
# in order to keep a fixed z_equality(R/M) and z_equality(M/Lambda)
#
omega_b = 0.0223828
omega_cdm_standard = 0.1201075
h_standard = 0.67810
#
# coefficient such that omega_r = omega_gamma (1 + coeff*Neff),
# i.e. such that omega_ur = omega_gamma * coeff * Neff:
# coeff = omega_ur/omega_gamma/Neff_standard 
# We could extract omega_ur and omega_gamma on-the-fly within th script, 
# but for simplicity we did a preliminary interactive run with background_verbose=2
# and we copied the values given in the budget output.
#
coeff = 1.70961e-05/2.47298e-05/3.044
print ("coeff=",coeff)
#
#############################################
#
# Fixed settings
#
############################################
#
# Varying parameter (others fixed to default)
#
var_name = 'N_ur'
var_array = [3.044]
var_num = len(var_array)
var_legend = r'$N_\mathrm{eff}$'
var_figname = 'neff'
#
# Constraints to be matched
#
# As explained in the "Neutrino cosmology" book, CUP, Lesgourgues et al., section 5.3, the goal is to vary
# - omega_cdm by a factor alpha = (1 + coeff*Neff)/(1 + coeff*3.046)
# - h by a factor sqrt*(alpha)
# in order to keep a fixed z_equality(R/M) and z_equality(M/Lambda)
#
omega_b = 0.0223828
omega_cdm_standard = 0.1201075
h_standard = 0.67810
#
# coefficient such that omega_r = omega_gamma (1 + coeff*Neff),
# i.e. such that omega_ur = omega_gamma * coeff * Neff:
# coeff = omega_ur/omega_gamma/Neff_standard 
# We could extract omega_ur and omega_gamma on-the-fly within th script, 
# but for simplicity we did a preliminary interactive run with background_verbose=2
# and we copied the values given in the budget output.
#
coeff = 1.70961e-05/2.47298e-05/3.044
print ("coeff=",coeff)
#
#############################################
#
# Fixed settings
#

common_settings = {'A_s':2.101e-9,
          'n_s':0.9665,
          'tau_reio':0.0561,
          'omega_b':0.02242,
          'omega_cdm':0.11933,
          'h':0.6766,
          'YHe':0.2425,
          'T_cmb':2.7255,
          'gauge':'newtonian', #FOR MGCLASS TO WORK, GAUGE NEEDS TO BE NEWTONIAN
          'k_pivot': 0.05,
          'mg_z_init': 10.000,
          'l_logstep': 1.025,
          'l_linstep':15,
          'P_k_max_1/Mpc':113.0,
          'l_switch_limber':9,
          'perturb_sampling_stepsize': 0.05,
          'output':'tCl,pCl,lCl,mPk',
          'l_max_scalars': 3000,
          'lensing': 'yes',
          'mg_ansatz':'z_xpans_late'}

#
##############################################
#
# loop over varying parameter values
#
beta_arr = [-1,0,1]
K0_arr = np.linspace(-0.75,1,10)


labels = [r'$g_{\mu}=0$   $g_{\gamma}=0$',
          r'$g_{\mu}=0.1$ $g_{\gamma}=0.1$',
          r'$g_{\mu}=0.1$ $g_{\gamma}=0.2$',
          r'$g_{\mu}=0.1$ $g_{\gamma}=0.5$',
          r'$g_{\mu}=0.2$ $g_{\gamma}=0.1$',
          r'$g_{\mu}=0.2$ $g_{\gamma}=0.2$',
          r'$g_{\mu}=0.2$ $g_{\gamma}=0.5$',
          r'$g_{\mu}=0.5$ $g_{\gamma}=0.1$',
          r'$g_{\mu}=0.5$ $g_{\gamma}=0.2$',
          r'$g_{\mu}=0.5$ $g_{\gamma}=0.5$']

M = {}
#
for j, beta in enumerate(beta_arr):
    for i, K0 in enumerate(K0_arr):
        common_settings['mg_T1'] = beta
        common_settings['mg_T2'] = K0
        common_settings['mg_T3'] = beta
        common_settings['mg_T4'] = K0
        common_settings['mg_zzn'] = 1

        #
        # rescale omega_cdm and h
        #
        #
        # call CLASS
        #
        M = Class()
        M.set(common_settings)

        M.compute()
        
        


        kvec = np.logspace(-4,np.log10(113),1000) # array of kvec in h/Mpc
        twopi = 2.*math.pi
        #
        # Create figures
        #

        #
        # loop over varying parameter values
        #
        ll = {}
        clM = {}
        clTT = {}
        pkM = {}
        legarray = []
        Tcmb = 2.72e6

        #
        #
        # deal with colors and legends
        #
        h = 0.6766
        #
        # get Cls
        #
        # store P(k) for common k values
        #
        pkM = []
        # The function .pk(k,z) wants k in 1/Mpc so we must convert kvec for each case with the right h 
        khvec = kvec*h # This is k in 1/Mpc
        for kh in khvec:
            pkM.append(M.pk(kh,0.)*h**3) 
        #    
        # plot P(k)
        #
        f = mticker.ScalarFormatter(useMathText=True)
        f.set_powerlimits((-6,6))
        ax_Pk.plot(kvec,pkM,
                        color=colors[j][i],#alpha=var_alpha,
                        linestyle='-')

    #
    # plot C_l^TT
    #
    
x = np.loadtxt('wmap.txt')[:,0]
y = np.loadtxt('wmap.txt')[:,1]
z = np.loadtxt('wmap.txt')[:,2]
plines_full = []
color = 'k'
plines =ax_Pk.errorbar(x,y,yerr=(z-y)*2,capsize=2,ecolor=color,color='w',marker='^',markersize=4,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 3)

plines_full.append(plines)

x = np.loadtxt('lya.txt')[:,0]
y = np.loadtxt('lya.txt')[:,1]
z = np.loadtxt('lya.txt')[:,2]

color = 'tab:gray'
plines = ax_Pk.errorbar(x,y,yerr=(z-y)*2,capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 3)
plines_full.append(plines)

legend1 = ax_Pk.legend(plines_full, [r"$\rm  WMAP/ACT$",  r"$\rm \mathrm{Ly-}\alpha\rm \,forest$"], loc='lower left',fancybox=True, fontsize=11)
legend1.get_frame().set_facecolor('none')
legend1.get_frame().set_linewidth(0.0)
ax_Pk.add_artist(legend1)

ax_Pk.set_xscale('log')
ax_Pk.set_yscale('log')
#plt.ylim(10**1.6,10**4.75)
#plt.xlim(10**(-3),1)
#plt.legend(loc='best')
ax_Pk.grid(".")
"""
h, l = ax_Pk.get_legend_handles_labels()
kw = dict(ncol=2, loc="lower center",fancybox=True, fontsize=11,frameon=False)    
leg2 = ax_Pk.legend(h[:],l[:], bbox_to_anchor=[0.5,1.08],**kw)
"""
ax_Pk.set_xscale('log')
ax_Pk.set_yscale('log')

ax_Pk.set_xlim(10**(-4),4*10**0)
ax_Pk.set_ylim(10**(0),5*10**5)


ax_Pk.set_xlabel(r'$k\;[h\,\mathrm{Mpc}^{-1}]$', size = '16')
ax_Pk.set_ylabel(r'$P_{\rm m}(k)\;[(h^{-1}\mathrm{Mpc})^3]$', size = '16')



#plt.legend(loc='best')
plt.grid(".")
line3 = ax_Pk.plot([0], [0], label=r'$T_{1}=-1$', color='#66c2a5')   
line2 = ax_Pk.plot([0], [0], label=r'$T_{1}=0$', color='#fc8d62')
line1 = ax_Pk.plot([0], [0], label=r'$T_{1}=1$', color='#8da0cb')

legend1 = ax_Pk.legend(loc='upper right',fancybox=True, fontsize=10)
legend1.get_frame().set_facecolor('none')
legend1.get_frame().set_linewidth(0.0)
ax_Pk.add_artist(legend1)


pl.rcParams['font.family'] = 'sans-serif'
p2 = ax_Pk.get_position().get_points().flatten()
#norm = colorss.Norm(pars1.min(), pars1.max())
K0_arr = np.array(K0_arr)
norm = mpl.colors.Normalize(vmin=K0_arr.min(), vmax=K0_arr.max())
x = 0.0167
ax_cbar = fig.add_axes([p2[0]-0.0218, 0.425, p2[2]-p2[0]+0.085, 0.01])
cbar_ax = plt.colorbar(mpl.cm.ScalarMappable(cmap=mpl.cm.Greys, norm=norm), cax=ax_cbar, orientation='horizontal', location = 'bottom', ticks=LinearLocator(numticks=8))
cbar_ax.set_label(r'$T_2$', fontsize=16)
cbar_ax.ax.tick_params(width=1.5, length=5, which = 'major')
cbar_ax.ax.tick_params(width=1.1, length=4, which = 'minor')
cbar_ax.ax.xaxis.set_minor_locator(AutoMinorLocator())
cbar_ax.ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(r'$'+'%.2g'+r'$'))
labels = cbar_ax.ax.get_xticklabels()
labels[0] = ""
cbar_ax.ax.set_xticklabels(labels)

ax_Pk = plt.subplot(424)
n = 10
cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#66c2a5"]) 
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#fc8d62"]) 
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#8da0cb"])

cmap = np.array([cmap1, cmap2, cmap3])
colors = [None]*3
colors[0] = cmap[0]((np.linspace(0, 1, n)))
colors[1] = cmap[1]((np.linspace(0, 1, n)))
colors[2] = cmap[2]((np.linspace(0, 1, n)))


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


############################################
#
# Varying parameter (others fixed to default)
#
var_name = 'N_ur'
var_array = [3.044]
var_num = len(var_array)
var_legend = r'$N_\mathrm{eff}$'
var_figname = 'neff'
#
# Constraints to be matched
#
# As explained in the "Neutrino cosmology" book, CUP, Lesgourgues et al., section 5.3, the goal is to vary
# - omega_cdm by a factor alpha = (1 + coeff*Neff)/(1 + coeff*3.046)
# - h by a factor sqrt*(alpha)
# in order to keep a fixed z_equality(R/M) and z_equality(M/Lambda)
#
omega_b = 0.0223828
omega_cdm_standard = 0.1201075
h_standard = 0.67810
#
# coefficient such that omega_r = omega_gamma (1 + coeff*Neff),
# i.e. such that omega_ur = omega_gamma * coeff * Neff:
# coeff = omega_ur/omega_gamma/Neff_standard 
# We could extract omega_ur and omega_gamma on-the-fly within th script, 
# but for simplicity we did a preliminary interactive run with background_verbose=2
# and we copied the values given in the budget output.
#
coeff = 1.70961e-05/2.47298e-05/3.044
print ("coeff=",coeff)
#
#############################################
#
# Fixed settings
#
############################################
#
# Varying parameter (others fixed to default)
#
var_name = 'N_ur'
var_array = [3.044]
var_num = len(var_array)
var_legend = r'$N_\mathrm{eff}$'
var_figname = 'neff'
#
# Constraints to be matched
#
# As explained in the "Neutrino cosmology" book, CUP, Lesgourgues et al., section 5.3, the goal is to vary
# - omega_cdm by a factor alpha = (1 + coeff*Neff)/(1 + coeff*3.046)
# - h by a factor sqrt*(alpha)
# in order to keep a fixed z_equality(R/M) and z_equality(M/Lambda)
#
omega_b = 0.0223828
omega_cdm_standard = 0.1201075
h_standard = 0.67810
#
# coefficient such that omega_r = omega_gamma (1 + coeff*Neff),
# i.e. such that omega_ur = omega_gamma * coeff * Neff:
# coeff = omega_ur/omega_gamma/Neff_standard 
# We could extract omega_ur and omega_gamma on-the-fly within th script, 
# but for simplicity we did a preliminary interactive run with background_verbose=2
# and we copied the values given in the budget output.
#
coeff = 1.70961e-05/2.47298e-05/3.044
print ("coeff=",coeff)
#
#############################################
#
# Fixed settings
#

common_settings = {'A_s':2.101e-9,
          'n_s':0.9665,
          'tau_reio':0.0561,
          'omega_b':0.02242,
          'omega_cdm':0.11933,
          'h':0.6766,
          'YHe':0.2425,
          'T_cmb':2.7255,
          'gauge':'newtonian', #FOR MGCLASS TO WORK, GAUGE NEEDS TO BE NEWTONIAN
          'k_pivot': 0.05,
          'mg_z_init': 10.000,
          'l_logstep': 1.025,
          'l_linstep':15,
          'P_k_max_1/Mpc':113.0,
          'l_switch_limber':9,
          'perturb_sampling_stepsize': 0.05,
          'output':'tCl,pCl,lCl,mPk',
          'l_max_scalars': 3000,
          'lensing': 'yes',
          'mg_ansatz':'GI',
          'Omega_Lambda': 0.0}


#
##############################################
#
# loop over varying parameter values
#
beta_arr = [-1.5,-1,-0.5]
K0_arr = np.linspace(0.4,0.6,10)


labels = [r'$g_{\mu}=0$   $g_{\gamma}=0$',
          r'$g_{\mu}=0.1$ $g_{\gamma}=0.1$',
          r'$g_{\mu}=0.1$ $g_{\gamma}=0.2$',
          r'$g_{\mu}=0.1$ $g_{\gamma}=0.5$',
          r'$g_{\mu}=0.2$ $g_{\gamma}=0.1$',
          r'$g_{\mu}=0.2$ $g_{\gamma}=0.2$',
          r'$g_{\mu}=0.2$ $g_{\gamma}=0.5$',
          r'$g_{\mu}=0.5$ $g_{\gamma}=0.1$',
          r'$g_{\mu}=0.5$ $g_{\gamma}=0.2$',
          r'$g_{\mu}=0.5$ $g_{\gamma}=0.5$']

M = {}
#
for j, beta in enumerate(beta_arr):
    for i, K0 in enumerate(K0_arr):
        common_settings['w0_fld'] = beta
        common_settings['gamGI'] = K0


        #
        # rescale omega_cdm and h
        #
        #
        # call CLASS
        #
        M = Class()
        M.set(common_settings)

        M.compute()
        
        


        kvec = np.logspace(-4,np.log10(113),1000) # array of kvec in h/Mpc
        twopi = 2.*math.pi
        #
        # Create figures
        #

        #
        # loop over varying parameter values
        #
        ll = {}
        clM = {}
        clTT = {}
        pkM = {}
        legarray = []
        Tcmb = 2.72e6

        #
        #
        # deal with colors and legends
        #
        h = 0.6766
        #
        # get Cls
        #
        # store P(k) for common k values
        #
        pkM = []
        # The function .pk(k,z) wants k in 1/Mpc so we must convert kvec for each case with the right h 
        khvec = kvec*h # This is k in 1/Mpc
        for kh in khvec:
            pkM.append(M.pk(kh,0.)*h**3) 
        #    
        # plot P(k)
        #
        f = mticker.ScalarFormatter(useMathText=True)
        f.set_powerlimits((-6,6))
        ax_Pk.plot(kvec,pkM,
                        color=colors[j][i],#alpha=var_alpha,
                        linestyle='-')

    #
    # plot C_l^TT
    #
    
x = np.loadtxt('wmap.txt')[:,0]
y = np.loadtxt('wmap.txt')[:,1]
z = np.loadtxt('wmap.txt')[:,2]
plines_full = []
color = 'k'
plines =ax_Pk.errorbar(x,y,yerr=(z-y)*2,capsize=2,ecolor=color,color='w',marker='^',markersize=4,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 3)

plines_full.append(plines)

x = np.loadtxt('lya.txt')[:,0]
y = np.loadtxt('lya.txt')[:,1]
z = np.loadtxt('lya.txt')[:,2]

color = 'tab:gray'
plines = ax_Pk.errorbar(x,y,yerr=(z-y)*2,capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 3)
plines_full.append(plines)

legend1 = ax_Pk.legend(plines_full, [r"$\rm  WMAP/ACT$",  r"$\rm \mathrm{Ly-}\alpha\rm \,forest$"], loc='lower left',fancybox=True, fontsize=11)
legend1.get_frame().set_facecolor('none')
legend1.get_frame().set_linewidth(0.0)
ax_Pk.add_artist(legend1)

ax_Pk.set_xscale('log')
ax_Pk.set_yscale('log')
#plt.ylim(10**1.6,10**4.75)
#plt.xlim(10**(-3),1)
#plt.legend(loc='best')
ax_Pk.grid(".")
"""
h, l = ax_Pk.get_legend_handles_labels()
kw = dict(ncol=2, loc="lower center",fancybox=True, fontsize=11,frameon=False)    
leg2 = ax_Pk.legend(h[:],l[:], bbox_to_anchor=[0.5,1.08],**kw)
"""
ax_Pk.set_xscale('log')
ax_Pk.set_yscale('log')

ax_Pk.set_xlim(10**(-4),4*10**0)
ax_Pk.set_ylim(10**(0),5*10**5)


ax_Pk.set_xlabel(r'$k\;[h\,\mathrm{Mpc}^{-1}]$', size = '16')



#plt.legend(loc='best')
plt.grid(".")
line3 = ax_Pk.plot([0], [0], label=r'$w_{\Lambda}=-1.5$', color='#66c2a5')   
line2 = ax_Pk.plot([0], [0], label=r'$w_{\Lambda}=-1$', color='#fc8d62')
line1 = ax_Pk.plot([0], [0], label=r'$w_{\Lambda}=-0.5$', color='#8da0cb')

legend1 = ax_Pk.legend(loc='upper right',fancybox=True, fontsize=10)
legend1.get_frame().set_facecolor('none')
legend1.get_frame().set_linewidth(0.0)
ax_Pk.add_artist(legend1)


pl.rcParams['font.family'] = 'sans-serif'
p2 = ax_Pk.get_position().get_points().flatten()
#norm = colorss.Norm(pars1.min(), pars1.max())
K0_arr = np.array(K0_arr)
norm = mpl.colors.Normalize(vmin=K0_arr.min(), vmax=K0_arr.max())
x = 0.0167
ax_cbar = fig.add_axes([p2[0]-0.0218+x, 0.425, p2[2]-p2[0]+0.085, 0.01])
cbar_ax = plt.colorbar(mpl.cm.ScalarMappable(cmap=mpl.cm.Greys, norm=norm), cax=ax_cbar, orientation='horizontal', location = 'bottom', ticks=LinearLocator(numticks=8))
cbar_ax.set_label(r'$\gamma$', fontsize=16)
cbar_ax.ax.tick_params(width=1.5, length=5, which = 'major')
cbar_ax.ax.tick_params(width=1.1, length=4, which = 'minor')
cbar_ax.ax.xaxis.set_minor_locator(AutoMinorLocator())
cbar_ax.ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(r'$'+'%.2g'+r'$'))
labels = cbar_ax.ax.get_xticklabels()
labels[0] = ""
cbar_ax.ax.set_xticklabels(labels)

ax_Pk.yaxis.set_ticklabels([])

plt.tight_layout()

plt.subplots_adjust(wspace=0, hspace=0)

plt.savefig('Pk_All.pdf',bbox_inches='tight')