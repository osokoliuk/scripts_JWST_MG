from JWST_MG.cosmological_functions import cosmological_library
model = 'nDGP'
model_H = 'nDGP'

par1 = 15000
par2 = 1

library = cosmological_library(1, model, model_H, par1, par2)
z = 2
print(library.dH_f(1/(1+z), model_H, par1, par2))

