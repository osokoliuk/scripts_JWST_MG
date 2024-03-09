from JWST_MG.UVLF import UVLF
from JWST_MG.HMF import HMF
from JWST_MG.reionization import reionization
from JWST_MG.delta_c import delta_c

from JWST_MG.constants import *

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
cbar = plt.colorbar(mpl.cm.ScalarMappable(cmap=pl.cm.Blues, norm=norm), ax=ax)
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
plt.savefig('Delta_vir_pheno.pdf', bbox_inches='tight')
