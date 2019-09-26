import auxiliary_functions as af
from ReadOutputGaussian import ReadOutputGaussian
import pandas as pd
import numpy as np



n_en = 35 
mut_in = "/home/mana/Desktop/python/test_ciro_pot/ci_mut_gamess_35.inp"
pot_in = "/home/mana/Desktop/python/test_ciro_pot/ci_pot_gamess_35.inp"

mut_out = "/home/mana/Desktop/python/test_ciro_pot/ci_mut_35.inp"
pot_out = "/home/mana/Desktop/python/test_ciro_pot/ci_pot_35.inp"




read = ReadOutputGaussian()

#mut
mut_gamess = pd.read_csv(mut_in, header=None, sep=r"\s+")
mut_val = np.array(mut_gamess.iloc[:,4:])
mut_name = np.array(mut_gamess.iloc[:,:4])

il=np.tril_indices(n_en)
il_index = np.array(il)

zero = np.hstack((mut_name[:n_en+1], mut_val[:n_en+1]))
reordered_val = read.read_half_below_matrix_gaussian(n_en+1,3,mut_val)
reordered_val = reordered_val[1:,1:]
reordered_val_tril = reordered_val[il]
reordered_names = list()
for (i,j) in zip(il_index[0,:],il_index[1,:]):
    reordered_names.append("States " + str(i+1) + " and " + str(j+1))

reordered_names = np.array(reordered_names)
mut_nozero = np.vstack((reordered_names.T, reordered_val_tril[:,0], reordered_val_tril[:,1], reordered_val_tril[:,2])).T

f = open(mut_out, 'w')
f.write("#ordered transition dipoles \n")
f.close()
f = open(mut_out, 'ab')
np.savetxt(f, zero, delimiter=' ', fmt="%s")
np.savetxt(f,mut_nozero, delimiter=' ', fmt="%s")
f.close()


pot_gamess = read.convertVgamess_to_VWaveT(pot_in,pot_out,n_en+1)







