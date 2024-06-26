

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

def compute_cumulative_uv_luminosity(muv_arr, phi_uv_arr, mlim_lum = -29):
    """
    Compute the cumulative UV luminosity
    assuming the Muv grid is uniform
    """
    dmuv = np.abs(muv_arr[1] - muv_arr[0]) # a constant
    Luv = 10**((muv_arr + 48.6)/(-2.5)) * 4*np.pi* (10*con.pc.value*100)**2  # erg/s/Hz
    cum_rho_uv_arr = np.zeros(len(phi_uv_arr))
    for i in range(len(phi_uv_arr)):
        select = (muv_arr[i:] > mlim_lum)
        cum_rho_uv_arr[i] = np.sum( (phi_uv_arr[i:] * Luv[i:])[select]  ) * dmuv
    return cum_rho_uv_arr # [#erg/s/Hz/Mpc^3]

def compute_total_uv_luminosity(muv_arr, phi_uv_arr, mlim = -17):
    """
    Compute the total UV luminosity
    assuming the Muv grid is uniform
    """
    dmuv = np.abs(muv_arr[1] - muv_arr[0]) # a constant
    select = (muv_arr < mlim) & (muv_arr > -29)
    Luv = 10**((muv_arr + 48.6)/(-2.5)) * 4*np.pi* (10*con.pc.value*100)**2  # erg/s/Hz
    return np.sum(phi_uv_arr[select] * Luv[select]) * dmuv # [#erg/s/Hz/Mpc^3]

def compute_cumulative_uv_count(muv_arr, phi_uv_arr):
    """
    Compute the cumulative UV count
    assuming the Muv grid is uniform
    """
    dmuv = np.abs(muv_arr[1] - muv_arr[0]) # a constant
    cum_phi_uv_arr = np.zeros(len(phi_uv_arr))
    for i in range(len(phi_uv_arr)):
        select = (muv_arr[i:] > -28)
        cum_phi_uv_arr[i] = np.sum(phi_uv_arr[i:][select]) * dmuv
    return cum_phi_uv_arr # [#/Mpc^3]

def compute_total_uv_count(muv_arr, phi_uv_arr):
    """
    Compute the total UV luminosity
    assuming the Muv grid is uniform
    """
    dmuv = np.abs(muv_arr[1] - muv_arr[0]) # a constant
    select = (muv_arr < -12) & (muv_arr > -28)
    return np.sum(phi_uv_arr[select]) * dmuv # [#/Mpc^3]

def compute_mean_shift(sigma, sigma_ref=1):
    #sigma is in unit of mag
    # shift in mean of a log-normal distribution
    sigmaL      = sigma/2.5
    sigmaL_ref  = sigma_ref/2.5
    meanX     = np.exp(( sigmaL * np.log(10) )**2 / 2.)
    meanX_ref = np.exp(( sigmaL_ref * np.log(10) )**2 / 2.)
    return 2.5 * np.log10(meanX/meanX_ref)

#print(compute_mean_shift(np.array([1,2,3,4])))

def load_harikane2023_specz(redshift, type=0):
    """
    Load the data from Harikane et al. 2023 
    """
    f = np.genfromtxt('../observational_data/Harikane2023_Specz.dat', names=True, delimiter=',')
    select = (f['z'] == redshift) & (f['type'] == type)
    log_density = np.log10(f['density'][select])
    density_lower = np.array([max(1e-10, x) for x in (f['density'][select] - f['lower_err'][select])])

    log_lower_err = log_density - np.log10(density_lower)
    log_upper_err = np.log10(f['density'][select] + f['upper_err'][select]) - log_density

    return f['ABmag'][select], log_density, log_lower_err, log_upper_err

def plot_specz_constraints(redshift, ax=None, **kwargs):
    x,y,ylo,yup = load_harikane2023_specz(redshift=redshift, type=0)
    print(y,ylo,yup)
    ax.errorbar( x, 10**y, yerr=(10**y-10**(y-ylo), -10**y+10**(y+yup)), **kwargs)
    x,y,ylo,yup = load_harikane2023_specz(redshift=redshift, type=1)
    ax.errorbar( x, 10**y, yerr=(10**y-10**(y-ylo), -10**y+10**(y+yup)), lolims=True, **kwargs)
    ax.errorbar( x, 10**y, yerr=(10**y-10**(y-ylo), -10**y+10**(y+yup)), **kwargs)
    x,y,ylo,yup = load_harikane2023_specz(redshift=redshift, type=2)
    ylo = np.array([0.3 for t in ylo if t==0])
    ax.errorbar( x, 10**y, yerr=(10**y-10**(y-ylo), -10**y+10**(y+yup)), uplims=True,**kwargs)
    ax.errorbar( x, 10**y, yerr=(10**y-10**(y-ylo), -10**y+10**(y+yup)), **kwargs)


def load_harikane_2023_photoz_digit(redshift):
    """
    Load the data from Harikane et al. 2023 
    """
    f = np.genfromtxt('../observational_data/Harikane2023_Photoz_digit.dat', names=True)
    select = (f['z'] == redshift)
    log_phi_err = (f["logPhi"][select] - f["lo"][select],  f["hi"][select] - f["logPhi"][select])
    mag_err     = (f["ABmag"][select] - f["left"][select],  f["right"][select] - f["ABmag"][select])

    return f['ABmag'][select], f["logPhi"][select], mag_err, log_phi_err

def plot_photoz_constraints(redshift, ax=None, **kwargs):
    basepath = "./observational_data/UVLF_photoz/"

    ########################################################
    f = np.genfromtxt(basepath + "Bouwens2021+Oesch2018.dat" , names=True)
    sel = (f['z'] == redshift) & (f['Phi'] - f['dPhi'] > 0)
    x = f['Muv'][sel]; y = np.log10(f['Phi'][sel])  
    yup = np.log10(f['Phi'][sel] + f['dPhi'][sel])  
    ylo = np.log10(f['Phi'][sel] - f['dPhi'][sel])
    ax.errorbar(x, 10**y, yerr=(10**(y-ylo), 10**(yup-y)), **kwargs)

    sel = (f['z'] == redshift) & (f['Phi'] - f['dPhi'] == 0)
    x = f['Muv'][sel]; y = np.log10(f['Phi'][sel])  
    yup = np.log10(f['Phi'][sel] + f['dPhi'][sel]) 
    ax.errorbar(x, 10**y, yerr=10**(0.3 * np.ones(len(y)), yup-y), uplims=True, **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Bowler2020.dat" , names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]; y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ylo = np.log10(f['Phi'][sel] - f['lo'][sel])
    ax.errorbar(x, 10**y, xerr=f['dMuv'][sel], yerr=(10**(y-ylo), 10**(yup-y)), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Castellano2023.dat" , names=True)
    sel = (f['z'] == redshift) & (f['type'] == 0)
    x = f['Muv'][sel]; y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ylo = np.log10(f['Phi'][sel] - f['lo'][sel])
    ax.errorbar(x, 10**y, xerr=f['dMuv'][sel], yerr=(10**(y-ylo), 10**(yup-y)), **kwargs)

    sel = (f['z'] == redshift) & (f['type'] == 2)
    x = f['Muv'][sel]; y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ax.errorbar(x, 10**y, xerr=f['dMuv'][sel], yerr=10**(0.3*np.ones(len(y)), yup-y), uplims=True, **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Donnan2023.dat" , names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]; y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ylo = np.log10(f['Phi'][sel] - f['lo'][sel])
    ax.errorbar(x, 10**y, xerr=f['dMuv2'][sel]/2., yerr=(10**(y-ylo), 10**(yup-y)), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Finkelstein2022a.dat" , names=True)
    if redshift == 10:
        y = np.log10(f['Phi'])
        yup = np.log10(f['Phi'] + f['up']); ylo = np.log10(f['Phi'] - f['lo'])
        ax.fill_between(10**f['Muv'], 10**ylo, 10**yup, color='gray', edgecolor='gray', alpha=0.5)

    ########################################################
    f = np.genfromtxt(basepath + "Finkelstein2022b.dat" , names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]; y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ylo = np.log10(f['Phi'][sel] - f['lo'][sel])
    ax.errorbar(x, 10**y, yerr=(10**(y-ylo), 10**(yup-y)), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Harikane2022.dat" , names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]; y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ylo = np.log10(f['Phi'][sel] - f['lo'][sel])
    ax.errorbar(x, 10**y, xerr=f['dMuv'][sel], yerr=(10**(y-ylo), 10**(yup-y)), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Harikane2023.dat" , names=True)
    sel = (f['z'] == redshift) & (f['type'] == 0)
    x = f['Muv'][sel]; y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ylo = np.log10(f['Phi'][sel] - f['lo'][sel])
    ax.errorbar(x, 10**y, xerr=f['dMuv'][sel], yerr=(10**(y-ylo), 10**(yup-y)), **kwargs)

    sel = (f['z'] == redshift) & (f['type'] == 2)
    x = f['Muv'][sel]; y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ax.errorbar(x, y, xerr=f['dMuv'][sel], yerr=(0.3*np.ones(len(y)), yup-y), uplims=True, **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "McLeod2016.dat" , names=True)
    sel = f['z'] == redshift
    x=f['Muv'][sel]; y=f['logPhi'][sel]
    yup = f['up'][sel]; ylo = f['lo'][sel] 
    xl  = f['left'][sel]; xr = f['right'][sel]
    ax.errorbar(x, y, xerr=(x-xl, xr-x), yerr=(10**(y-ylo), 10**(yup-y)), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Morishita2018+2023.dat" , names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]; y=f['logPhi'][sel]; dyup=f['up'][sel]; dylo=f['lo'][sel]
    ax.errorbar(x, y, yerr=(dylo, dyup), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Naidu2023.dat" , names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]; y=f['logPhi'][sel]; dyup=f['up'][sel]; dylo=f['lo'][sel]
    ax.errorbar(x, 10**y, xerr=f["dMuv"][sel], yerr=(10**dylo, 10**dyup), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Perez2023.dat" , names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]; y=f['logPhi'][sel]; dyup=f['up'][sel]; dylo=f['lo'][sel]
    ax.errorbar(x, 10**y, yerr=(10**dylo, 10**dyup), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Stefanon2019.dat" , names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]; y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ylo = np.log10(f['Phi'][sel] - f['lo'][sel])
    ax.errorbar(x, 10**y, yerr=(10**(y-ylo), 10**(yup-y)), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "McLeod2023.dat" , names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]; y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['dPhi'][sel])
    ylo = np.log10(f['Phi'][sel] - f['dPhi'][sel])
    ax.errorbar(x, 10**y, yerr=(10**(y-ylo), 10**(yup-y)), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Adams2023.dat" , names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]; y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['dPhi'][sel])
    ylo = np.log10(f['Phi'][sel] - f['dPhi'][sel])
    ax.errorbar(x, 10**y, yerr=(10**(y-ylo), 10**(yup-y)), **kwargs)
    

    
plt.cla()
plt.figure()
plt.rcParams.update({"text.usetex": True})
fig = plt.figure(figsize=(4.25*2*.95*1.05, 2*4*1.05*1.05))


nn = 1
z_smf_arr = [4]

for z_smf in z_smf_arr:
    ax_Pk = plt.subplot(4,2,1)

    
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



    obs = number_density(feature='GLF_UV', z_target=z_smf, h=h)
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
                ax_Pk.errorbar(data[:,0],  data[:,1],yerr = np.abs([data[:,1]-data[:,3],data[:,2]- data[:,1]]),\
                        label=r'$\rm pre-JWST$',capsize=0,ecolor=color,color='w',marker=marker,markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color, zorder= 3)
            else:
                ax_Pk.errorbar(data[:,0],  data[:,1],yerr = np.abs([data[:,1]-data[:,3],data[:,2]- data[:,1]]),\
                        capsize=0,ecolor=color,color='w',marker=marker,markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color, zorder= 3)
            
            j_data +=1


    """if z_smf in [4,5,6,7,8]:
        path = '../observational_data/GSMF'
        Navarro = np.loadtxt(path + "/Navarro_z"+str(z_smf)+".dat")
        x = 10**Navarro[:,0]
        y = 1e-4*Navarro[:,1]
        yerr = 1e-4*Navarro[:,2]
        color = 'k'
        ax_Pk.errorbar(x,y,yerr=yerr,label=r'$\rm JWST$',capsize=0,ecolor=color,color='w',marker='v',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color, zorder= 3)
    """


    color = 'k'
    marker = 'v'
    plot_specz_constraints(redshift=int(z_smf), ax=ax_Pk,capsize=0,ecolor=color,color='w',marker=marker,markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color, zorder= 3)

    if z_smf in [12,10,9]:
        ax_Pk.errorbar(-100, 1, yerr=1, capsize=0,ecolor=color,color='w',marker=marker,markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color, zorder= 3, label = r'$\rm JWST$')

    pool_cpu = Pool(8)

    model = 'nDGP'
    model_H = 'nDGP'
    model_SFR = 'double_power'
    pars1 = np.logspace(2.5, 5, 10)
    par2 = 0
    f0 = 0.21
    n = len(pars1)
    cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#48639e"])

    colors = cmap3(np.linspace(0, 1, n))    
    Pk_arr = []
    for par1 in pars1:
        HMF_library = HMF(1/(1+z_smf), model, model_H, par1, par2, 1e8)
        Pk_arr.append(np.array(HMF_library.Pk(1/(1+z_smf), model, par1, par2))*h**3)
    k = kvec/h
    Masses = np.logspace(5,19,250)

    UVLF_library = UVLF(1/(1+z_smf), model, model_H, model_SFR, pars1, par2, Masses, f0)

    sigma_uv = 0.4
    iterable = [(1/(1+z_smf), rhom, model, model_H, model_SFR, par1, par2, Masses, k, Pk_arr[i], f0, sigma_uv) for i,par1 in enumerate(pars1)]
    MUV, UVLF_obs = zip(*pool_cpu.starmap(UVLF_library.compute_uv_luminosity_function,tqdm(iterable, total=len(pars1))))
    for i in range(len(UVLF_obs)):
        ax_Pk.plot(MUV[i], UVLF_obs[i], c = colors[i], lw=  1)

    norm = colorss.LogNorm(pars1.min(), pars1.max())
    cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap3, norm=norm), ax=ax_Pk)
    cbar.set_label(r'$r_c$', fontsize=16)

    pars1 = np.array([1e8])
    Pk_arr = []
    for par1 in pars1:
        HMF_library = HMF(1/(1+z_smf), model, model_H, par1, par2, 1e8)
        Pk_arr.append(np.array(HMF_library.Pk(1/(1+z_smf), model, par1, par2))*h**3)
    k = kvec/h
    Masses = np.logspace(5,19,250)

    UVLF_library = UVLF(1/(1+z_smf), model, model_H, model_SFR, pars1, par2, Masses, f0)

    sigmas = np.array([2,4])
    iterable = [(1/(1+z_smf), rhom, model, model_H, model_SFR, par1, par2, Masses, k, Pk_arr[0], f0, sigma_uv) for sigma_uv in sigmas]
    MUV, UVLF_obs = zip(*pool_cpu.starmap(UVLF_library.compute_uv_luminosity_function,tqdm(iterable, total=len(sigmas))))
    for i in range(len(UVLF_obs)):
        line, = ax_Pk.plot(MUV[i], UVLF_obs[i], c = 'k', lw=4, alpha=0.2)
        line_annotate(r'$\sigma_{\rm UV}=' + str(sigmas[i]) + '$',line,-23, c = 'tab:gray', fontsize = 9)



    #plt.errorbar(x.get('Duncan'),y.get('Duncan'),yerr=[yerr_down.get('Duncan'),yerr_up.get('Duncan')], c = 'tab:orange', capsize = 2, ls = 'None', marker = '.', label = r'$\rm Duncan+14$')
    #plt.errorbar(x.get('Song'),y.get('Song'),yerr=[yerr_down.get('Song'),yerr_up.get('Song')], c = 'tab:orange', capsize = 2, ls = 'None', marker = 's', label = r'$\rm Song+16$')
    #plines = plt.errorbar(x.get('Duncan'),y.get('Duncan'),yerr=[yerr_down.get('Duncan'),yerr_up.get('Duncan')],capsize=0,ecolor='tab:blue',color='w',marker='o',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='tab:blue')
    #plines = plt.errorbar(x.get('Song'),y.get('Song'),yerr=[yerr_down.get('Song'),yerr_up.get('Song')],capsize=0,ecolor='tab:orange',color='w',marker='s',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='tab:orange')


    #plines = plt.errorbar(x.get('Navarro'),y.get('Navarro'),yerr=[yerr_down.get('Navarro'),yerr_up.get('Navarro')],capsize=0,ecolor='k',color='w',marker=markers[j_data+1],markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='k', label = r'$\rm Navarro+2024$')

    # plt.scatter(1/a_vir-1, vir2, c = 'tab:orange')
    #plt.xscale('log')
    plt.yscale('log')
    plt.xlim(-26,-12)
    plt.ylim(10**(-8),10**(-0.5))

    ax_Pk.set_xticklabels([])
    ax_Pk.set_ylabel(r'$\phi_{\rm UV}\;[\rm Mpc^{-3}\;mag^{-1}]$', size = '16')


    plt.grid(".")
    
    ax_Pk.text(-25,10**(-1.75),r'$z='+str(int(round(z_smf)))+r'$', size = '15')
    
    legend1 = ax_Pk.legend(loc='lower right',fancybox=True, fontsize=10)
    legend1.get_frame().set_facecolor('none')
    legend1.get_frame().set_linewidth(0.0)
    ax_Pk.add_artist(legend1)
    
    nn += 1


nn = 2
z_smf_arr = [5,6,7,8,9,10,12]

for z_smf in z_smf_arr:
    ax_Pk = plt.subplot(4,2,nn)

    
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


    obs = number_density(feature='GLF_UV', z_target=z_smf, h=h)
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
                ax_Pk.errorbar(data[:,0],  data[:,1],yerr = np.abs([data[:,1]-data[:,3],data[:,2]- data[:,1]]),\
                        label=r'$\rm pre-JWST$',capsize=0,ecolor=color,color='w',marker=marker,markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color, zorder= 3)
            else:
                ax_Pk.errorbar(data[:,0],  data[:,1],yerr = np.abs([data[:,1]-data[:,3],data[:,2]- data[:,1]]),\
                        capsize=0,ecolor=color,color='w',marker=marker,markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color, zorder= 3)
            
            j_data +=1


    color = 'k'
    marker = 'v'
    plot_specz_constraints(redshift=int(z_smf), ax=ax_Pk,capsize=0,ecolor=color,color='w',marker=marker,markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color, zorder= 3)

    if z_smf in [12,10,9]:
        ax_Pk.errorbar(-100, 1, yerr=1, capsize=0,ecolor=color,color='w',marker=marker,markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor=color, zorder= 3, label = r'$\rm JWST$')

    model = 'nDGP'
    model_H = 'nDGP'
    model_SFR = 'double_power'
    pars1 = np.logspace(2.5, 5, 10)
    par2 = 0
    f0 = 0.21
    n = len(pars1)
    cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white","#48639e"])

    colors = cmap3(np.linspace(0, 1, n))    
    Pk_arr = []
    for par1 in pars1:
        HMF_library = HMF(1/(1+z_smf), model, model_H, par1, par2, 1e8)
        Pk_arr.append(np.array(HMF_library.Pk(1/(1+z_smf), model, par1, par2))*h**3)
    k = kvec/h
    Masses = np.logspace(5,19,250)

    UVLF_library = UVLF(1/(1+z_smf), model, model_H, model_SFR, pars1, par2, Masses, f0)

    sigma_uv = 0.4
    iterable = [(1/(1+z_smf), rhom, model, model_H, model_SFR, par1, par2, Masses, k, Pk_arr[i], f0, sigma_uv) for i,par1 in enumerate(pars1)]
    MUV, UVLF_obs = zip(*pool_cpu.starmap(UVLF_library.compute_uv_luminosity_function,tqdm(iterable, total=len(pars1))))
    for i in range(len(UVLF_obs)):
        ax_Pk.plot(MUV[i], UVLF_obs[i], c = colors[i], lw=  1)

    norm = colorss.LogNorm(pars1.min(), pars1.max())
    cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap3, norm=norm), ax=ax_Pk)
    cbar.set_label(r'$r_c$', fontsize=16)

    pars1 = np.array([1e8])
    Pk_arr = []
    for par1 in pars1:
        HMF_library = HMF(1/(1+z_smf), model, model_H, par1, par2, 1e8)
        Pk_arr.append(np.array(HMF_library.Pk(1/(1+z_smf), model, par1, par2))*h**3)
    k = kvec/h
    Masses = np.logspace(5,19,250)

    UVLF_library = UVLF(1/(1+z_smf), model, model_H, model_SFR, pars1, par2, Masses, f0)

    sigmas = np.array([2,4])
    iterable = [(1/(1+z_smf), rhom, model, model_H, model_SFR, par1, par2, Masses, k, Pk_arr[0], f0, sigma_uv) for sigma_uv in sigmas]
    MUV, UVLF_obs = zip(*pool_cpu.starmap(UVLF_library.compute_uv_luminosity_function,tqdm(iterable, total=len(sigmas))))
    for i in range(len(UVLF_obs)):
        line, = ax_Pk.plot(MUV[i], UVLF_obs[i], c = 'k', lw=4, alpha=0.2)
        if nn == 8:
            line_annotate(r'$\sigma_{\rm UV}=' + str(sigmas[i]) + '$',line,-22.5, c = 'tab:gray', fontsize = 9)
        else:
            line_annotate(r'$\sigma_{\rm UV}=' + str(sigmas[i]) + '$',line,-23, c = 'tab:gray', fontsize = 9)
        

    #plt.errorbar(x.get('Duncan'),y.get('Duncan'),yerr=[yerr_down.get('Duncan'),yerr_up.get('Duncan')], c = 'tab:orange', capsize = 2, ls = 'None', marker = '.', label = r'$\rm Duncan+14$')
    #plt.errorbar(x.get('Song'),y.get('Song'),yerr=[yerr_down.get('Song'),yerr_up.get('Song')], c = 'tab:orange', capsize = 2, ls = 'None', marker = 's', label = r'$\rm Song+16$')
    #plines = plt.errorbar(x.get('Duncan'),y.get('Duncan'),yerr=[yerr_down.get('Duncan'),yerr_up.get('Duncan')],capsize=0,ecolor='tab:blue',color='w',marker='o',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='tab:blue')
    #plines = plt.errorbar(x.get('Song'),y.get('Song'),yerr=[yerr_down.get('Song'),yerr_up.get('Song')],capsize=0,ecolor='tab:orange',color='w',marker='s',markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='tab:orange')


    #plines = plt.errorbar(x.get('Navarro'),y.get('Navarro'),yerr=[yerr_down.get('Navarro'),yerr_up.get('Navarro')],capsize=0,ecolor='k',color='w',marker=markers[j_data+1],markersize=4,markeredgewidth=1, elinewidth=1.2,ls='None',markeredgecolor='k', label = r'$\rm Navarro+2024$')

    # plt.scatter(1/a_vir-1, vir2, c = 'tab:orange')
    #plt.xscale('log')
    plt.yscale('log')
    plt.xlim(-26,-12)
    plt.ylim(10**(-8),10**(-0.5))
    if nn == 7 or nn == 8:
        ax_Pk.set_xlabel(r'$M_{\rm UV}\;[\rm mag]$', size = '16')
        ax_Pk.set_ylabel(r'$\phi_{\rm UV}\;[\rm Mpc^{-3}\;mag^{-1}]$', size = '16')
    else:
        ax_Pk.set_xticklabels([])
        ax_Pk.set_ylabel(r'$\phi_{\rm UV}\;[\rm Mpc^{-3}\;mag^{-1}]$', size = '16')


    plt.grid(".")
    
    ax_Pk.text(-25,10**(-1.75),r'$z='+str(int(round(z_smf)))+r'$', size = '15')
    
    legend1 = ax_Pk.legend(loc='lower right',fancybox=True, fontsize=10)
    legend1.get_frame().set_facecolor('none')
    legend1.get_frame().set_linewidth(0.0)
    ax_Pk.add_artist(legend1)
    
    nn += 1


plt.tight_layout()

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
3
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


plt.tight_layout()
plt.savefig('UVLF_nDGP_double_power.pdf', bbox_inches='tight')
