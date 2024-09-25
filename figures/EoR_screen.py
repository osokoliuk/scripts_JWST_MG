

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
from matplotlib.text import Annotation
from matplotlib.transforms import Affine2D

class LineAnnotation(Annotation):
    """A sloped annotation to *line* at position *x* with *text*
    Optionally an arrow pointing from the text to the graph at *x* can be drawn.
    Usage
    -----
    fig, ax = subplots()
    x = linspace(0, 2*pi)
    line, = ax.plot(x, sin(x))
    ax.add_artist(LineAnnotation("text", line, 1.5))
    """

    def __init__(
        self, text, line, x, xytext=(0, 5), textcoords="offset points", **kwargs
    ):
        """Annotate the point at *x* of the graph *line* with text *text*.

        By default, the text is displayed with the same rotation as the slope of the
        graph at a relative position *xytext* above it (perpendicularly above).

        An arrow pointing from the text to the annotated point *xy* can
        be added by defining *arrowprops*.

        Parameters
        ----------
        text : str
            The text of the annotation.
        line : Line2D
            Matplotlib line object to annotate
        x : float
            The point *x* to annotate. y is calculated from the points on the line.
        xytext : (float, float), default: (0, 5)
            The position *(x, y)* relative to the point *x* on the *line* to place the
            text at. The coordinate system is determined by *textcoords*.
        **kwargs
            Additional keyword arguments are passed on to `Annotation`.

        See also
        --------
        `Annotation`
        `line_annotate`
        """
        assert textcoords.startswith(
            "offset "
        ), "*textcoords* must be 'offset points' or 'offset pixels'"

        self.line = line
        self.xytext = xytext

        # Determine points of line immediately to the left and right of x
        xs, ys = line.get_data()

        def neighbours(x, xs, ys, try_invert=True):
            inds, = np.where((xs <= x)[:-1] & (xs > x)[1:])
            if len(inds) == 0:
                assert try_invert, "line must cross x"
                return neighbours(x, xs[::-1], ys[::-1], try_invert=False)

            i = inds[0]
            return np.asarray([(xs[i], ys[i]), (xs[i+1], ys[i+1])])
        
        self.neighbours = n1, n2 = neighbours(x, xs, ys)
        
        # Calculate y by interpolating neighbouring points
        y = n1[1] + ((x - n1[0]) * (n2[1] - n1[1]) / (n2[0] - n1[0]))

        kwargs = {
            "horizontalalignment": "center",
            "rotation_mode": "anchor",
            **kwargs,
        }
        super().__init__(text, (x, y), xytext=xytext, textcoords=textcoords, **kwargs)

    def get_rotation(self):
        """Determines angle of the slope of the neighbours in display coordinate system
        """
        transData = self.line.get_transform()
        dx, dy = np.diff(transData.transform(self.neighbours), axis=0).squeeze()
        return np.rad2deg(np.arctan2(dy, dx))

    def update_positions(self, renderer):
        """Updates relative position of annotation text
        Note
        ----
        Called during annotation `draw` call
        """
        xytext = Affine2D().rotate_deg(self.get_rotation()).transform(self.xytext)
        self.set_position(xytext)
        super().update_positions(renderer)


def line_annotate(text, line, x, *args, **kwargs):
    """Add a sloped annotation to *line* at position *x* with *text*

    Optionally an arrow pointing from the text to the graph at *x* can be drawn.

    Usage
    -----
    x = linspace(0, 2*pi)
    line, = ax.plot(x, sin(x))
    line_annotate("sin(x)", line, 1.5)

    See also
    --------
    `LineAnnotation`
    `plt.annotate`
    """
    ax = line.axes
    a = LineAnnotation(text, line, x, *args, **kwargs)
    if "clip_on" in kwargs:
        a.set_clip_path(ax.patch)
    ax.add_artist(a)
    return a


    
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
fig = plt.figure(figsize=(4.25*2*.95*0.9, 2*5*1.05*0.9))


nn = 1
z_smf_arr = [4,5,6,7,8,9]
turn_tau = False

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


    ax_Pk.spines['top'].set_zorder(5)
    ax_Pk.spines['right'].set_zorder(5)
    ax_Pk.spines['left'].set_zorder(5)
    ax_Pk.spines['bottom'].set_zorder(5)


    if nn % 2 != 0:
        data = np.loadtxt('../observational_data/EoR/QHII.txt')
        x = data[:,0]
        y = data[:,1]
        yerr = data[:,2]
        lims_switch = data[:,3]
        marker = '.'
        color = 'tab:gray'
        plt.errorbar(x,1-y,yerr=yerr, lolims=lims_switch,capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 4)


        if nn < 3:
            model = 'nDGP'
            model_H = 'nDGP'
            model_SFR = 'Puebla'
            pars1 = np.logspace(2.5,6,10)
            #pars1 = np.hstack((pars1,[10**4.5,10**5]))
            par2 = 0
            f0 = 0.21
            n = len(pars1)
            cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#fc8d62"])

            colors = cmap3(np.linspace(0, 1, n))
            
            
            for i, par1 in enumerate(pars1):
                colors = cmap3(np.linspace(0, 1, n))
                #ac_arr = np.linspace(1/21, 1, 64)
                #reion = reionization(ac_arr, model, model_H, par1, par2)
                #z_int, qhii = reion.QHII(rhom, model, model_H, model_SFR, par1, par2, pool_cpu, f0=f0)
                #data = np.savez('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+'.npz', name1 = z_int, name2 = qhii)
                data = np.load('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+'.npz')
                z_int = data['name1']
                qhii = data['name2']
                #print(reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0))
                #tau = reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0)
                #print(tau)
                ax_Pk.plot(z_int,qhii, c = colors[i], lw=  1.5, zorder = 3)
                #pool_cpu.terminate()


            model = 'nDGP'
            model_H = 'nDGP'
            model_SFR = 'double_power'
            pars1 = np.logspace(2.5,6,10)
            #pars1 = np.hstack((pars1,[10**4.5,10**5]))
            par2 = 0
            f0 = 0.21
            n = len(pars1)
            cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#48639e"])

            colors = cmap3(np.linspace(0, 1, n))
            
            
            for i, par1 in enumerate(pars1):
                colors = cmap3(np.linspace(0, 1, n))
                #ac_arr = np.linspace(1/21, 1, 64)
                #reion = reionization(ac_arr, model, model_H, par1, par2)
                #z_int, qhii = reion.QHII(rhom, model, model_H, model_SFR, par1, par2, pool_cpu, f0=f0)
                #data = np.savez('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+'.npz', name1 = z_int, name2 = qhii)
                data = np.load('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+'.npz')
                z_int = data['name1']
                qhii = data['name2']
                #print(reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0))
                #tau = reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0)
                #print(tau)
                ax_Pk.plot(z_int,qhii, c = colors[i], lw=  1.5, zorder = 3)
                #pool_cpu.terminate()
                
                
            ax_Pk.plot(0,0,c = '#fc8d62', label = r'$\rm Rodriguez-Puebla$')
            ax_Pk.plot(0,0,c = '#8da0cb', label = r'$\rm Double\;power-law$')
            # plt.scatter(1/a_vir-1, vir2, c = 'tab:orange')
            plt.xlim(4,12)
            plt.ylim(0,1)
        if nn == 3:
            data = np.loadtxt('../observational_data/EoR/QHII.txt')
            x = data[:,0]
            y = data[:,1]
            yerr = data[:,2]
            lims_switch = data[:,3]
            marker = '.'
            color = 'tab:gray'
            #plt.errorbar(x,1-y,yerr=yerr, lolims=lims_switch,capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 3)
            
            model = 'kmoufl'
            model_H = 'kmoufl'
            model_SFR = 'Puebla'
            pars2 = np.linspace(0,1,10) #np.array([0.1])
            pars1 = np.array([0.1, 0.3, 0.5]) #np.array([0.1])
            f0 = 0.21
            n = len(pars2)

            cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#66c2a5"]) 
            cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#fc8d62"]) 
            cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#8da0cb"])

            cmap = np.array([cmap1, cmap2, cmap3])
            colors = [None]*3
            colors[0] = cmap[0]((np.linspace(0, 1, n)))
            colors[1] = cmap[1]((np.linspace(0, 1, n)))
            colors[2] = cmap[2]((np.linspace(0, 1, n)))
                        
            for i, par1 in enumerate(pars1):
                for j, par2 in enumerate(pars2):
                    #ac_arr = np.linspace(1/21, 1, 64)
                    #reion = reionization(ac_arr, model, model_H, par1, par2)
                    #z_int, qhii = reion.QHII(rhom, model, model_H, model_SFR, par1, par2, pool_cpu, f0=f0)
                    #data = np.savez('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+'.npz', name1 = z_int, name2 = qhii)
                    try:
                        data = np.load('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+ '_' +str(par2)+'.npz')
                        z_int = data['name1']
                        qhii = data['name2']
                        #print(reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0))
                        #tau = reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0)
                        #print(tau)
                        ax_Pk.plot(z_int,qhii, c = colors[i][j], lw=  1.5, zorder = 3)
                        #pool_cpu.terminate()
                    except:
                        print('failed')
            ax_Pk.plot(0,0, lw=  1.5, zorder = 3, c = colors[0][-1],label = r"$\beta = 0.1$")
            ax_Pk.plot(0,0, lw=  1.5, zorder = 3, c = colors[1][-1],label = r"$\beta = 0.3$")
            ax_Pk.plot(0,0, lw=  1.5, zorder = 3, c = colors[2][-1],label = r"$\beta = 0.5$")

            plt.xlim(4,12)
            plt.ylim(0,1)
        elif nn == 5:
            data = np.loadtxt('../observational_data/EoR/QHII.txt')
            x = data[:,0]
            y = data[:,1]
            yerr = data[:,2]
            lims_switch = data[:,3]
            marker = '.'
            color = 'tab:gray'
            #plt.errorbar(x,1-y,yerr=yerr, lolims=lims_switch,capsize=2,ecolor=color,color='w',marker=marker,markersize=6,markeredgewidth=1.3, elinewidth=1,ls='None',markeredgecolor=color, zorder= 3)
            
            model = 'kmoufl'
            model_H = 'kmoufl'
            model_SFR = 'double_power'
            pars2 = np.linspace(0,1,10) #np.array([0.1])
            pars1 = np.array([0.1, 0.3, 0.5]) #np.array([0.1])
            f0 = 0.21
            n = len(pars2)

            cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#66c2a5"]) 
            cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#fc8d62"]) 
            cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#8da0cb"])

            cmap = np.array([cmap1, cmap2, cmap3])
            colors = [None]*3
            colors[0] = cmap[0]((np.linspace(0, 1, n)))
            colors[1] = cmap[1]((np.linspace(0, 1, n)))
            colors[2] = cmap[2]((np.linspace(0, 1, n)))
                        
            for i, par1 in enumerate(pars1):
                for j, par2 in enumerate(pars2):
                    #ac_arr = np.linspace(1/21, 1, 64)
                    #reion = reionization(ac_arr, model, model_H, par1, par2)
                    #z_int, qhii = reion.QHII(rhom, model, model_H, model_SFR, par1, par2, pool_cpu, f0=f0)
                    #data = np.savez('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+'.npz', name1 = z_int, name2 = qhii)
                    try:
                        data = np.load('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+ '_' +str(par2)+'.npz')
                        z_int = data['name1']
                        qhii = data['name2']
                        #print(reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0))
                        #tau = reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0)
                        #print(tau)
                        ax_Pk.plot(z_int,qhii, c = colors[i][j], lw=  1.5, zorder = 3)
                        #pool_cpu.terminate()
                    except:
                        print('failed')
            ax_Pk.plot(0,0, lw=  1.5, zorder = 3, c = colors[0][-1],label = r"$\beta = 0.1$")
            ax_Pk.plot(0,0, lw=  1.5, zorder = 3, c = colors[1][-1],label = r"$\beta = 0.3$")
            ax_Pk.plot(0,0, lw=  1.5, zorder = 3, c = colors[2][-1],label = r"$\beta = 0.5$")

            plt.xlim(4,12)
            plt.ylim(0,1)
    else:
        if nn < 3:
            if turn_tau == True:
                model = 'nDGP'
                model_H = 'nDGP'
                model_SFR = 'Puebla'
                pars1 = np.logspace(2.5,6,10)
                #pars1 = np.hstack((pars1,[10**4.5,10**5]))
                par2 = 0
                f0 = 0.21
                n = len(pars1)
                cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#fc8d62"])

                colors = cmap3(np.linspace(0, 1, n))
                
                
                for i, par1 in enumerate(pars1):
                    colors = cmap3(np.linspace(0, 1, n))
                    ac_arr = np.linspace(1/21, 1, 64)
                    reion = reionization(ac_arr, model, model_H, par1, par2)
                    #z_int, qhii = reion.QHII(rhom, model, model_H, model_SFR, par1, par2, pool_cpu, f0=f0)
                    #data = np.savez('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+'.npz', name1 = z_int, name2 = qhii)
                    data = np.load('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+'.npz')
                    z_int = data['name1']
                    qhii = data['name2'][:,0]
                    qhii_int = scipy.interpolate.interp1d(z_int,qhii,fill_value='extrapolate')
                    tau = []
                    for j, z_end in enumerate(z_int):
                        z_span = np.linspace(z_end,0,100)
                        tau.append(reion.tau_reio(z_span, qhii_int(z_span), model, model_H, par1, par2))
                    #print(reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0))
                    #tau = reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0)
                    #print(tau)
                    ax_Pk.plot(z_int,tau, c = colors[i], lw=  1.5, zorder = 3)
                    #pool_cpu.terminate()
                
                model = 'nDGP'
                model_H = 'nDGP'
                model_SFR = 'double_power'
                pars1 = np.logspace(2.5,6,10)
                #pars1 = np.hstack((pars1,[10**4.5,10**5]))
                par2 = 0
                f0 = 0.21
                n = len(pars1)
                cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#48639e"])

                colors = cmap3(np.linspace(0, 1, n))
                
                
                for i, par1 in enumerate(pars1):
                    colors = cmap3(np.linspace(0, 1, n))
                    ac_arr = np.linspace(1/21, 1, 64)
                    reion = reionization(ac_arr, model, model_H, par1, par2)
                    #z_int, qhii = reion.QHII(rhom, model, model_H, model_SFR, par1, par2, pool_cpu, f0=f0)
                    #data = np.savez('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+'.npz', name1 = z_int, name2 = qhii)
                    data = np.load('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+'.npz')
                    z_int = data['name1']
                    qhii = data['name2'][:,0]
                    qhii_int = scipy.interpolate.interp1d(z_int,qhii,fill_value='extrapolate')
                    tau = []
                    for j, z_end in enumerate(z_int):
                        z_span = np.linspace(z_end,0,100)
                        tau.append(reion.tau_reio(z_span, qhii_int(z_span), model, model_H, par1, par2))
                    #print(reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0))
                    #tau = reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0)
                    #print(tau)
                    ax_Pk.plot(z_int,tau, c = colors[i], lw=  1.5, zorder = 3)
                    #pool_cpu.terminate()"""
            
                        
            plt.axhline(0.054, c='tab:gray', lw=0.8)
            plt.text(2, 0.054-0.007-0.005, r'$\rm Planck\; 2018$',
            fontsize=11, c='tab:grey')
            ax_Pk.fill_between([0,12], 0.054-0.007, 0.054 +
                    0.007, alpha=0.25, color='tab:gray')

            ax_Pk.plot(0,0,c = '#fc8d62', label = r'$\rm Rodriguez-Puebla$')
            ax_Pk.plot(0,0,c = '#8da0cb', label = r'$\rm Double\;power-law$')
            ax_Pk.yaxis.tick_right()
            ax_Pk.yaxis.set_label_position("right")
            plt.xlim(0,12)
            plt.ylim(0,0.07)
        elif nn == 4:
            if turn_tau == True:
                model = 'kmoufl'
                model_H = 'kmoufl'
                model_SFR = 'double_power'
                pars2 = np.linspace(0,1,10) #np.array([0.1])
                pars1 = np.array([0.1, 0.3, 0.5]) #np.array([0.1])
                f0 = 0.21
                n = len(pars2)

                cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#66c2a5"]) 
                cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#fc8d62"]) 
                cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#8da0cb"])

                cmap = np.array([cmap1, cmap2, cmap3])
                colors = [None]*3
                colors[0] = cmap[0]((np.linspace(0, 1, n)))
                colors[1] = cmap[1]((np.linspace(0, 1, n)))
                colors[2] = cmap[2]((np.linspace(0, 1, n)))
                            
                for i, par1 in enumerate(pars1):
                    for j, par2 in enumerate(pars2):
                        ac_arr = np.linspace(1/21, 1, 120)
                        reion = reionization(ac_arr, model, model_H, par1, par2)
                        #z_int, qhii = reion.QHII(rhom, model, model_H, model_SFR, par1, par2, pool_cpu, f0=f0)
                        #data = np.savez('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+'.npz', name1 = z_int, name2 = qhii)
                        data = np.load('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+ '_' +str(par2)+'.npz')
                        z_int = data['name1']
                        qhii = data['name2'][:,0]
                        qhii_int = scipy.interpolate.interp1d(z_int,qhii,fill_value='extrapolate')
                        tau = []
                        for l, z_end in enumerate(z_int):
                            z_span = np.linspace(z_end,0,300)
                            tau.append(reion.tau_reio(z_span, qhii_int(z_span), model, model_H, par1, par2))
                        #print(reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0))
                        #tau = reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0)
                        #print(tau)
                        ax_Pk.plot(z_int,tau, c = colors[i][j], lw=  1.5, zorder = 3)
                        #pool_cpu.terminate()"""
            
                        
            plt.axhline(0.054, c='tab:gray', lw=0.8)
            plt.text(2, 0.054-0.007-0.005, r'$\rm Planck\; 2018$',
            fontsize=11, c='tab:grey')
            ax_Pk.fill_between([0,12], 0.054-0.007, 0.054 +
                    0.007, alpha=0.25, color='tab:gray')

            ax_Pk.plot(0,0, lw=  1.5, zorder = 3, c = colors[0][-1],label = r"$\beta = 0.1$")
            ax_Pk.plot(0,0, lw=  1.5, zorder = 3, c = colors[1][-1],label = r"$\beta = 0.3$")
            ax_Pk.plot(0,0, lw=  1.5, zorder = 3, c = colors[2][-1],label = r"$\beta = 0.5$")

            ax_Pk.yaxis.tick_right()
            ax_Pk.yaxis.set_label_position("right")
            plt.xlim(0,12)
            plt.ylim(0,0.07)
        elif nn == 6:
            if turn_tau == True:
                model = 'kmoufl'
                model_H = 'kmoufl'
                model_SFR = 'double_power'
                pars2 = np.linspace(0,1,10) #np.array([0.1])
                pars1 = np.array([0.1, 0.3, 0.5]) #np.array([0.1])
                f0 = 0.21
                n = len(pars2)

                cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#66c2a5"]) 
                cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#fc8d62"]) 
                cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#8da0cb"])

                cmap = np.array([cmap1, cmap2, cmap3])
                colors = [None]*3
                colors[0] = cmap[0]((np.linspace(0, 1, n)))
                colors[1] = cmap[1]((np.linspace(0, 1, n)))
                colors[2] = cmap[2]((np.linspace(0, 1, n)))
                            
                for i, par1 in enumerate(pars1):
                    for j, par2 in enumerate(pars2):
                        ac_arr = np.linspace(1/21, 1, 120)
                        reion = reionization(ac_arr, model, model_H, par1, par2)
                        #z_int, qhii = reion.QHII(rhom, model, model_H, model_SFR, par1, par2, pool_cpu, f0=f0)
                        #data = np.savez('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+'.npz', name1 = z_int, name2 = qhii)
                        data = np.load('./data_folder/QHII_'+str(model) +'_'+ str(model_SFR) + '_' +str(par1)+ '_' +str(par2)+'.npz')
                        z_int = data['name1']
                        qhii = data['name2'][:,0]
                        qhii_int = scipy.interpolate.interp1d(z_int,qhii,fill_value='extrapolate')
                        tau = []
                        for l, z_end in enumerate(z_int):
                            z_span = np.linspace(z_end,0,300)
                            tau.append(reion.tau_reio(z_span, qhii_int(z_span), model, model_H, par1, par2))
                        #print(reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0))
                        #tau = reion.tau_reio(rhom, model, model_H, model_SFR, par1, par2, f0=f0)
                        #print(tau)
                        ax_Pk.plot(z_int,tau, c = colors[i][j], lw=  1.5, zorder = 3)
                        #pool_cpu.terminate()"""
            
                        
            plt.axhline(0.054, c='tab:gray', lw=0.8)
            plt.text(2, 0.054-0.007-0.005, r'$\rm Planck\; 2018$',
            fontsize=11, c='tab:grey')
            ax_Pk.fill_between([0,12], 0.054-0.007, 0.054 +
                    0.007, alpha=0.25, color='tab:gray')

            ax_Pk.plot(0,0, lw=  1.5, zorder = 3, c = colors[0][-1],label = r"$\beta = 0.1$")
            ax_Pk.plot(0,0, lw=  1.5, zorder = 3, c = colors[1][-1],label = r"$\beta = 0.3$")
            ax_Pk.plot(0,0, lw=  1.5, zorder = 3, c = colors[2][-1],label = r"$\beta = 0.5$")

            ax_Pk.yaxis.tick_right()
            ax_Pk.yaxis.set_label_position("right")
            plt.xlim(0,12)
            plt.ylim(0,0.07)
    if nn != 5 and nn != 6:
        if nn % 2 == 0:
            nbins = len(ax_Pk.get_yticklabels())
            ax_Pk.yaxis.set_major_locator(MaxNLocator(nbins=nbins,prune='lower'))
            ax_Pk.set_xticklabels([])
            #ax_Pk.set_yticklabels([])
            ax_Pk.set_ylabel(r'$\tau_{\rm reion}(z)$', size = '16')

        else:
            if nn == 1:
                nbins = len(ax_Pk.get_yticklabels())
                ax_Pk.yaxis.set_major_locator(MaxNLocator(nbins=nbins,prune='lower'))
            else:
                nbins = len(ax_Pk.get_yticklabels())
                ax_Pk.yaxis.set_major_locator(MaxNLocator(nbins=nbins,prune='both'))
            ax_Pk.set_xticklabels([])
            ax_Pk.set_ylabel(r'$Q_{\rm HII}(z)$', size = '16')
    else:
        if nn % 2 == 0:
            ax_Pk.set_xlabel(r'$\mathrm{Redshift}\;z$', size = '16')
            nbins = len(ax_Pk.get_yticklabels())
            ax_Pk.yaxis.set_major_locator(MaxNLocator(nbins=nbins,prune='upper'))
            nbins = len(ax_Pk.get_yticklabels())
            ax_Pk.yaxis.set_major_locator(MaxNLocator(nbins=nbins,prune='lower'))
            #ax_Pk.set_yticklabels([])
            ax_Pk.set_ylabel(r'$\tau_{\rm reion}(z)$', size = '16')

        else:
            ax_Pk.set_xlabel(r'$\mathrm{Redshift}\;z$', size = '16')
            ax_Pk.set_ylabel(r'$Q_{\rm HII}(z)$', size = '16')
            nbins = len(ax_Pk.get_yticklabels())
            ax_Pk.yaxis.set_major_locator(MaxNLocator(nbins=nbins,prune='upper'))
            ax_Pk.xaxis.set_major_locator(MaxNLocator(nbins=nbins,prune='upper'))
    #plt.grid(".")
    
    
    if nn % 2 == 0:
        legend1 = ax_Pk.legend(loc='lower right',fancybox=True, fontsize=10)
    else:
        legend1 = ax_Pk.legend(loc='upper right',fancybox=True, fontsize=10)
    legend1.get_frame().set_facecolor('none')
    legend1.get_frame().set_linewidth(0.0)
    ax_Pk.add_artist(legend1)
    
    nn += 1




mpl.rcParams['font.family'] = 'sans-serif'

#norm = colorss.Norm(pars1.min(), pars1.max())
norm = mpl.colors.Normalize(vmin=2.5, vmax=4)
ax_cbar = fig.add_axes([0.1275, 0.985, 0.8175, 0.01])
cbar_ax = plt.colorbar(mpl.cm.ScalarMappable(cmap=mpl.cm.Greys, norm=norm), cax=ax_cbar, orientation='horizontal', location = 'top', ticks=LinearLocator(numticks=8))
cbar_ax.set_label(r'$\log_{10}r_c$', fontsize=16)
cbar_ax.ax.tick_params(width=1.5, length=5, which = 'major')
cbar_ax.ax.tick_params(width=1.1, length=4, which = 'minor')
cbar_ax.ax.xaxis.set_minor_locator(AutoMinorLocator())
cbar_ax.ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(r'$'+'%.2g'+r'$'))


plt.tight_layout()

plt.subplots_adjust(wspace=0, hspace=0)

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
par2 = 0
ac_arr = np.linspace(0.01, 1, 15)

for ac in ac_arr:
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


#plt.tight_layout()
plt.savefig('EoR_screen.pdf', bbox_inches='tight')
