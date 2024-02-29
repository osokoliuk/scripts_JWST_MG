
from JWST_MG.UVLF import UVLF
from JWST_MG.HMF import HMF
from JWST_MG.reionization import reionization
from JWST_MG.delta_c import delta_c

from JWST_MG.constants import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import commah
ax = plt.subplot(111)
model = 'nDGP'
model_H = 'nDGP'
model_SFR = 'toy'

par1 = 150000
par2 = 1

Mh0 = 1e12
f0 = 0.05

######################################################################
plt.cla()

z = 10
Masses = np.logspace(8, 16, 50)
a = 1/(1+z)
ac = 1
reion = reionization(1, model, model_H, par1, par2)
deltac_library = delta_c(
    a_arr, model, model_H, par1, par2)
deltai = deltac_library.interpolate_ac(ac, model, model_H, par1, par2)
print(deltai)
delta_nl = deltac_library.non_linear(
    deltai, a_arr, model, model_H, par1, par2)
a_arr, R_arr = reion.radius_solve(model, model_H, par1, par2, deltai, delta_nl)
plt.plot(a_arr, R_arr)
plt.xlim(ai, ac)
plt.ylim(0, 7.5e4)
"""
def halo_accretion_rate(mhalo, redshift):
    # Fakhouri 2010
    # Mhalo in Msun
    # mhalo_dot = 46.1 * (1 + 1.11*redshift) * np.sqrt(Omegam0*(1+redshift)**3 + (1-Omegam0))  \
    # * (mhalo / 1e12)**(1.1)
    mhalo_dot = 25.3 * (1 + 1.65*redshift) * np.sqrt(Omegam0*(1+redshift)**3 + (1-Omegam0))  \
        * (mhalo / 1e12)**(1.1)
    # corr = 10**(-0.1) # down 0.1 dex consider the drop of sigma8
    return mhalo_dot


# dMhdt_mean = 25.3*(Mh0/1e12)**1.1*(1+1.65*z)*np.sqrt(Omegam0*(1+z)**3+1-Omegam0)
# plt.plot(np.log10(1+z), np.log10(halo_accretion_rate(Mh0, z)))


def load_harikane2023_specz(redshift, type=0):
    #Load the data from Harikane et al. 2023 
    f = np.genfromtxt(
        'observational_data/Harikane2023_Specz.dat', names=True, delimiter=',')
    select = (f['z'] == redshift) & (f['type'] == type)
    log_density = np.log10(f['density'][select])
    density_lower = np.array(
        [max(1e-10, x) for x in (f['density'][select] - f['lower_err'][select])])

    log_lower_err = log_density - np.log10(density_lower)
    log_upper_err = np.log10(f['density'][select] +
                             f['upper_err'][select]) - log_density

    return f['ABmag'][select], log_density, log_lower_err, log_upper_err


def plot_specz_constraints(redshift, ax=None, **kwargs):
    x, y, ylo, yup = load_harikane2023_specz(redshift=redshift, type=0)
    ax.errorbar(x, y, yerr=(ylo, yup), **kwargs)
    x, y, ylo, yup = load_harikane2023_specz(redshift=redshift, type=1)
    ax.errorbar(x, y, yerr=(ylo, yup), lolims=True, **kwargs)
    ax.errorbar(x, y, yerr=(ylo, yup), **kwargs)
    x, y, ylo, yup = load_harikane2023_specz(redshift=redshift, type=2)
    ylo = np.array([0.3 for t in ylo if t == 0])
    ax.errorbar(x, y, yerr=(ylo, yup), uplims=True, **kwargs)
    ax.errorbar(x, y, yerr=(ylo, yup), **kwargs)


def load_harikane_2023_photoz_digit(redshift):
    #Load the data from Harikane et al. 2023 
    f = np.genfromtxt(
        'observational_data/Harikane2023_Photoz_digit.dat', names=True)
    select = (f['z'] == redshift)
    log_phi_err = (f["logPhi"][select] - f["lo"][select],
                   f["hi"][select] - f["logPhi"][select])
    mag_err = (f["ABmag"][select] - f["left"][select],
               f["right"][select] - f["ABmag"][select])

    return f['ABmag'][select], f["logPhi"][select], mag_err, log_phi_err


def plot_photoz_constraints(redshift, ax=None, **kwargs):
    basepath = "./observational_data/UVLF_photoz/"

    ########################################################
    f = np.genfromtxt(basepath + "Bouwens2021+Oesch2018.dat", names=True)
    sel = (f['z'] == redshift) & (f['Phi'] - f['dPhi'] > 0)
    x = f['Muv'][sel]
    y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['dPhi'][sel])
    ylo = np.log10(f['Phi'][sel] - f['dPhi'][sel])
    ax.errorbar(x, y, yerr=(y-ylo, yup-y), **kwargs)

    sel = (f['z'] == redshift) & (f['Phi'] - f['dPhi'] == 0)
    x = f['Muv'][sel]
    y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['dPhi'][sel])
    ax.errorbar(x, y, yerr=(0.3 * np.ones(len(y)), yup-y),
                uplims=True, **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Bowler2020.dat", names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]
    y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ylo = np.log10(f['Phi'][sel] - f['lo'][sel])
    ax.errorbar(x, y, xerr=f['dMuv'][sel], yerr=(y-ylo, yup-y), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Castellano2023.dat", names=True)
    sel = (f['z'] == redshift) & (f['type'] == 0)
    x = f['Muv'][sel]
    y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ylo = np.log10(f['Phi'][sel] - f['lo'][sel])
    ax.errorbar(x, y, xerr=f['dMuv'][sel], yerr=(y-ylo, yup-y), **kwargs)

    sel = (f['z'] == redshift) & (f['type'] == 2)
    x = f['Muv'][sel]
    y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ax.errorbar(x, y, xerr=f['dMuv'][sel], yerr=(
        0.3*np.ones(len(y)), yup-y), uplims=True, **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Donnan2023.dat", names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]
    y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ylo = np.log10(f['Phi'][sel] - f['lo'][sel])
    ax.errorbar(x, y, xerr=f['dMuv2'][sel]/2., yerr=(y-ylo, yup-y), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Finkelstein2022a.dat", names=True)
    if redshift == 10:
        y = np.log10(f['Phi'])
        yup = np.log10(f['Phi'] + f['up'])
        ylo = np.log10(f['Phi'] - f['lo'])
        ax.fill_between(f['Muv'], ylo, yup, color='gray',
                        edgecolor='gray', alpha=0.5)

    ########################################################
    f = np.genfromtxt(basepath + "Finkelstein2022b.dat", names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]
    y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ylo = np.log10(f['Phi'][sel] - f['lo'][sel])
    ax.errorbar(x, y, yerr=(y-ylo, yup-y), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Harikane2022.dat", names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]
    y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ylo = np.log10(f['Phi'][sel] - f['lo'][sel])
    ax.errorbar(x, y, xerr=f['dMuv'][sel], yerr=(y-ylo, yup-y), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Harikane2023.dat", names=True)
    sel = (f['z'] == redshift) & (f['type'] == 0)
    x = f['Muv'][sel]
    y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ylo = np.log10(f['Phi'][sel] - f['lo'][sel])
    ax.errorbar(x, y, xerr=f['dMuv'][sel], yerr=(y-ylo, yup-y), **kwargs)

    sel = (f['z'] == redshift) & (f['type'] == 2)
    x = f['Muv'][sel]
    y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ax.errorbar(x, y, xerr=f['dMuv'][sel], yerr=(
        0.3*np.ones(len(y)), yup-y), uplims=True, **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "McLeod2016.dat", names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]
    y = f['logPhi'][sel]
    yup = f['up'][sel]
    ylo = f['lo'][sel]
    xl = f['left'][sel]
    xr = f['right'][sel]
    ax.errorbar(x, y, xerr=(x-xl, xr-x), yerr=(y-ylo, yup-y), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Morishita2018+2023.dat", names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]
    y = f['logPhi'][sel]
    dyup = f['up'][sel]
    dylo = f['lo'][sel]
    ax.errorbar(x, y, yerr=(dylo, dyup), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Naidu2023.dat", names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]
    y = f['logPhi'][sel]
    dyup = f['up'][sel]
    dylo = f['lo'][sel]
    ax.errorbar(x, y, xerr=f["dMuv"][sel], yerr=(dylo, dyup), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Perez2023.dat", names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]
    y = f['logPhi'][sel]
    dyup = f['up'][sel]
    dylo = f['lo'][sel]
    ax.errorbar(x, y, yerr=(dylo, dyup), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Stefanon2019.dat", names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]
    y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['up'][sel])
    ylo = np.log10(f['Phi'][sel] - f['lo'][sel])
    ax.errorbar(x, y, yerr=(y-ylo, yup-y), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "McLeod2023.dat", names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]
    y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['dPhi'][sel])
    ylo = np.log10(f['Phi'][sel] - f['dPhi'][sel])
    ax.errorbar(x, y, yerr=(y-ylo, yup-y), **kwargs)

    ########################################################
    f = np.genfromtxt(basepath + "Adams2023.dat", names=True)
    sel = f['z'] == redshift
    x = f['Muv'][sel]
    y = np.log10(f['Phi'][sel])
    yup = np.log10(f['Phi'][sel] + f['dPhi'][sel])
    ylo = np.log10(f['Phi'][sel] - f['dPhi'][sel])
    ax.errorbar(x, y, yerr=(y-ylo, yup-y), **kwargs)


# just for comparison at z<=8
def mv2020_plt_obdata(fname, papername, color, label=True):
    obdataPath = "observational_data/obdata_MV2020_archived/"
    obdata = np.genfromtxt(obdataPath+fname, names=True, comments='#')
    x_ob = obdata['m']
    phi = obdata['phi']
    id1 = phi > 0
    id2 = (phi <= 0)
    y_ob = 0.0*x_ob
    y_ob[id1] = np.log10(phi[id1])
    y_ob[id2] = phi[id2]
    uperr = 0.0*x_ob
    uperr[id1] = np.log10(phi[id1]+obdata['uperr'][id1])-y_ob[id1]
    uperr[id2] = obdata['uperr'][id2]
    low = phi-obdata['lowerr']
    low[low <= 0] = 1e-10
    lowerr = 0.0*x_ob
    lowerr[id1] = -np.log10(low[id1])+y_ob[id1]
    lowerr[id2] = obdata['lowerr'][id2]
    if label == True:
        ax.errorbar(x_ob, y_ob, yerr=(lowerr, uperr), c=color, lw=3,
                    linestyle='', marker='o', markersize=9, capsize=4.5, label=papername)
    else:
        ax.errorbar(x_ob, y_ob, yerr=(lowerr, uperr), c=color, lw=3,
                    linestyle='', marker='o', markersize=9, capsize=4.5)


redshift = 10
if redshift == 7:
    mv2020_plt_obdata('obdata'+str(11).zfill(3) +
                      '.dat', '', 'grey', label=False)
    mv2020_plt_obdata('obdata'+str(11).zfill(3)+'_ouc.dat',
                      r'${\rm Ouchi+}$ ${\rm 2009}$', 'saddlebrown')
    mv2020_plt_obdata('obdata'+str(11).zfill(3)+'_ate.dat',
                      r'${\rm Atek+}$ ${\rm 2015}$', 'silver')
    # plt_obdata('obdata'+str(Snaps[i]).zfill(3)+'_ono.dat','','lightgrey',label=False)
elif redshift == 6:
    mv2020_plt_obdata('obdata'+str(13).zfill(3) +
                      '.dat', '', 'grey', label=False)
    mv2020_plt_obdata('obdata'+str(13).zfill(3)+'_bou.dat',
                      r'${\rm Bouwens+}$ ${\rm 2017}$', 'lightgrey')
    mv2020_plt_obdata('obdata'+str(13).zfill(3)+'_ate.dat',
                      r'${\rm Atek+}$ ${\rm 2018}$', 'saddlebrown')
    # plt_obdata('obdata'+str(Snaps[i]).zfill(3)+'_ono.dat','','lightgrey',label=False)
elif redshift == 8:
    mv2020_plt_obdata('obdata'+str(8).zfill(3)+'.dat', '', 'grey', label=False)
elif redshift == 9:
    mv2020_plt_obdata('obdata'+str(4).zfill(3)+'.dat', '', 'grey', label=False)


my_ylims = np.r_[-8, -2]
my_xlims = np.r_[-16, -24.5]

kwargs = {"ms": 10,  "color": "silver", "mec": "silver", "mfc": "silver", "marker": "o", "linestyle": "none",
          "capsize": 4, "lw": 2, "mew": 2, "capthick": 2, "alpha": 1, "zorder": -8}
plot_photoz_constraints(redshift=redshift, ax=ax, **kwargs)
plt.errorbar([], [], label=r"$\rm JWST\,photo\,z$" +
             '\n' + r"$\rm +\,\,pre\,JWST$", **kwargs)

kwargs = {"ms": 13,  "color": "navy", "mec": "navy", "mfc": "navy", "marker": "o", "linestyle": "none",
          "capsize": 4, "lw": 2, "capthick": 2, "alpha": 1, "zorder": -8}
plot_specz_constraints(redshift=redshift, ax=ax, **kwargs)
plt.errorbar([], [], label=r"$\rm JWST\,spec\,z$", **kwargs)

leg = plt.legend(fontsize=22, loc=3, frameon=True, framealpha=1, numpoints=1)
leg.get_frame().set_linewidth(0)
plt.xlabel(r'$M_{\rm UV}$', fontsize=25)
plt.ylabel(r"$\log{(\Phi\,[{\rm mag}^{-1}\, {\rm Mpc}^{-3}])}$", fontsize=25)



my_ylims = np.r_[-8, -1]
my_xlims = np.r_[-14, -24]
plt.xlim(*my_xlims)
plt.ylim(*my_ylims)

"""

plt.savefig('HMF.pdf')
# delta_c_at_ac(self, ac, model, model_H, par1, par2):
