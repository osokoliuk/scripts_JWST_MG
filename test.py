from JWST_MG.UVLF import UVLF
from JWST_MG.HMF import HMF
from JWST_MG.reionization import reionization
from JWST_MG.delta_c import delta_c

from JWST_MG.constants import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

model = 'nDGP'
model_H = 'nDGP'
model_SFR = 'toy'

par2 = 6/11
par1 = 50000
deltac = delta_c(1, model, model_H, par1, par2)
print(deltac.binary_search_di(1, model, model_H, par1, par2, 0, 99))
