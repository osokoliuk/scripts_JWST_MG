
from JWST_MG.delta_c import delta_c
model = 'nDGP'
model_H = 'nDGP'

par1 = 150000
par2 = 1

delta_c_at_a1 = delta_c(1, model, model_H, par1, par2)
z = 15
print(delta_c_at_a1.delta_c_at_ac(1/(1+z),model, model_H, par1, par2))

#delta_c_at_ac(self, ac, model, model_H, par1, par2):
