import SystemManager as ini
import time
from multiprocessing import freeze_support
import numpy as np



#folder = "/Users/castd/Desktop/OCpy_hybrid_implementation/test/oc/scipy/3levels/"
folder = "/Users/castd/Desktop/OCpy_hybrid_implementation/"

distance = [distance for distance in np.arange(0.14, 2.14, 0.1)]

for bond_length in distance:

    namefile = "input.dat" + str(bond_length)




#OC_system=ini.SystemManager()
#OC_system.init_system(folder, namefile)

#print(OC_system.mol.par.muT[2,7])
#print(OC_system.mol.par.muT[7,1])

    if __name__ == "__main__":
        freeze_support()
    OC_system=ini.SystemManager()
    OC_system.init_system(folder, namefile)
    #start=time.time()
    OC_system.oc.iterate()
#end=time.time()
#print("serial: " + str(end-start))


#read = ReadOutputGaussian()
#read.convertVgamess_to_VWaveT(folder + "ci_pot.inp", folder + "ci_pot_transf.inp", 11)
