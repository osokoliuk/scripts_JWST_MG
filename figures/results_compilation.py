
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
from JWST_MG.UVLF import UVLF
import matplotlib.colors as colorss
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import LinearLocator

from mpl_toolkits.axes_grid1 import make_axes_locatable
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


    
import observation_data as obs

def MUV_lim(z):
    z_arr = [4,5,6,7,8]
    MUV_arr = [-22.6,-23,-22.5,-22.75,-22]
    MUV_lim_interp = scipy.interpolate.interp1d(z_arr,MUV_arr, fill_value='extrapolate')
    return MUV_lim_interp(z)



mpl.rcParams['axes.linewidth'] = 1.5

plt.cla()
plt.figure()
plt.rcParams.update({"text.usetex": True})
fig = plt.figure(figsize=(4.25*2*.95*0.9, 2*3*1.05*0.9))


ax = plt.subplot(1, 3, 1)
ax_Pk = ax
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


ax_Pk.spines['top'].set_zorder(5)
ax_Pk.spines['right'].set_zorder(5)
ax_Pk.spines['left'].set_zorder(5)
ax_Pk.spines['bottom'].set_zorder(5)

x = [3, 2.7, 3.9, 3, 3.9, (2.7 + 3)/2]
y = [5, 4, 2, 3, 1, 0]
xlo = [0.0, 1.0, 0.0, 0.0, 0.0, (3- 2.7)/2]
xup = [1, 0, 1, 2, 1, (3- 2.7)/2]
color = '#fc8d62'
marker = 'o'
ax.errorbar(x,y,yerr=None,xerr=[xlo,xup],capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 4)


x = [3.5, 2.9, 3, 2.8, 3, 3.5]
y = np.array([5, 4, 2, 3, 1, 0])+0.2
xlo = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
xup = [1, 0, 1, 2, 1, 1]
color = '#48639e'
marker = 'o'
ax.errorbar(x,y,yerr=None,xerr=[xlo,xup],capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 4)



ax.set_xlim(2.5,4)
ax.xaxis.set_major_locator(plt.MaxNLocator(4))
ax.set_ylim(-0.49,5.5)
ax.yaxis.set_major_locator(plt.MaxNLocator(6))
ax.set_xlabel(r'$\log_{10}r_c$', fontsize=16)
labels = [item.get_text() for item in ax.get_yticklabels()]
labels[6] = r'$\rm SMF$'
labels[5] = r'$\rm SMD$'
labels[4] = r'$\rm SFRD$'
labels[3] = r'$\rm UVLF$'
labels[2] = r'$Q_{\rm HII}$'
labels[1] = r'$\tau_{\rm reion}$'

ax.set_yticklabels(labels)

plt.fill_betweenx([-1,6], 3.5, 4, color='#48639e', alpha=.1)
plt.fill_betweenx([-1,6], 3.9, 4, color='#fc8d62', alpha=.1)
ax.yaxis.set_tick_params(labelsize=14)
ax.xaxis.set_tick_params(labelsize=12)

ax = plt.subplot(1, 3, 2)
ax_Pk = ax
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


ax_Pk.spines['top'].set_zorder(5)
ax_Pk.spines['right'].set_zorder(5)
ax_Pk.spines['left'].set_zorder(5)
ax_Pk.spines['bottom'].set_zorder(5)

x = [0.1, 0.3, 0.45, 0.1, 0.3, 0.3, 0.3, 0.45, 0.3, 0.45]
y = [5, 4, 4, 3, 3, 2, 1, 1, 0, 0]
xlo = [1, 0, 0, 1, 0, 0, 0, 0, 0, 0]
xup = [0, 0, 1, 0, 0, 0, 0, 1, 0, 1]
color = '#fc8d62'

marker_arr = np.array(['o','o','s','o','s','o','o','s','o','s'])
for i in range(len(marker_arr)):
    marker = marker_arr[i]
    ax.errorbar(x[i],y[i],yerr=None,xerr=[[xlo[i]],[xup[i]]],capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 4)



x = [0.1,0.3,0.45,0.1,0.1,0.1,0.1,0.3]
y = np.array([5,4,4,3,2,1,0,0]) + 0.2
xlo = [1,0,0,1,1,1,1,0]
xup = [0,0,1,0,0,0,0,0]
color = '#48639e'

marker_arr = np.array(['o','o','s','o','o','o','o','s'])
for i in range(len(marker_arr)):
    marker = marker_arr[i]
    ax.errorbar(x[i],y[i],yerr=None,xerr=[[xlo[i]],[xup[i]]],capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 4)




ax.set_xlim(0,0.5)
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_ylim(-0.49,5.5)
ax.yaxis.set_major_locator(plt.MaxNLocator(6))
ax.set_xlabel(r'$\beta$', fontsize=16)

ax.set_yticklabels([])

labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = ''
labels[-1] = ''

ax.set_xticklabels(labels)
ax.xaxis.set_tick_params(labelsize=12)

plt.fill_betweenx([-1,6], 0.3,0.5, color='#fc8d62', alpha=.1)
plt.fill_betweenx([-1,6], 0.0, 0.1, color='#48639e', alpha=.1)

ax = plt.subplot(1, 3, 3)
ax_Pk = ax
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


ax_Pk.spines['top'].set_zorder(5)
ax_Pk.spines['right'].set_zorder(5)
ax_Pk.spines['left'].set_zorder(5)
ax_Pk.spines['bottom'].set_zorder(5)

x = [0.9,0.5,0.8,0.5,0.9,0.5,0.5,0.9,(0.1+0.4)/2,(0.9+0.5)/2]
y = [5,4,4,3,3,2,1,1,0,0]
xlo = [0,1,0,1,0,1,1,0,(0.4-0.1)/2,(0.9-0.5)/2]
xup = [1,0,1,0,1,0,0,1,(0.4-0.1)/2,(0.9-0.5)/2]
color = '#fc8d62'

marker_arr = np.array(['o','o','s','o','s','o','o','s','o','s'])
for i in range(len(marker_arr)):
    marker = marker_arr[i]
    ax.errorbar(x[i],y[i],yerr=None,xerr=[[xlo[i]],[xup[i]]],capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 4)


x = [0.9,0.7,0.9,0.9,0.9,0.9,0.9,(0.9+0.7)/2]
y = (np.array([5,4,4,3,2,1,0,0]) + 0.2)[0:len(x)]
xlo = [0,1,0,0,0,0,0,(0.9-0.7)/2]
xup = [1,0,1,1,1,1,1,(0.9-0.7)/2]
color = '#48639e'
marker_arr = np.array(['o','o','s','o','o','o','o','s'])
for i in range(len(marker_arr)):
    marker = marker_arr[i]
    ax.errorbar(x[i],y[i],yerr=None,xerr=[[xlo[i]],[xup[i]]],capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 4)


ax.set_xlim(0,1)
ax.xaxis.set_major_locator(plt.MaxNLocator(4))
ax.set_ylim(-0.49,5.5)
ax.yaxis.set_major_locator(plt.MaxNLocator(6))
ax.set_xlabel(r'$K_0$', fontsize=16)

ax.set_yticklabels([])

labels = [item.get_text() for item in ax.get_xticklabels()]
ax.set_xticklabels(labels)
from matplotlib.ticker import StrMethodFormatter
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # No decimal places
ax.xaxis.set_tick_params(labelsize=12)

plt.fill_betweenx([-1,6], 0.0,0.3, color='#fc8d62', alpha=.1)
plt.fill_betweenx([-1,6], 0.9, 1.0, color='#48639e', alpha=.1)


plt.tight_layout()

plt.subplots_adjust(wspace=0, hspace=0)

marker = '.'
color = '#fc8d62'

ax.errorbar(-1,-1,xerr=0,capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 4,label=r'$\rm Rodriguez-Puebla$')
color = '#48639e'

ax.errorbar(-1,-1,xerr=0,capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 4,label=r'$\rm Double \;power-law$')



legend1 = plt.legend(bbox_to_anchor =(0.35,1), loc='lower right', ncols= 2,fancybox=True, fontsize = 11)
legend1.get_frame().set_facecolor('none')
legend1.get_frame().set_linewidth(0.0)

#plt.tight_layout()
plt.savefig('compiled.pdf', bbox_inches='tight')
