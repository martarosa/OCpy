import SystemManager as ini
import time
from multiprocessing import freeze_support


#folder = "/home/mana/programmi/python/optimal_control/OCpy/test_rabitz/"
#folder = "/home/mana/programmi/python/optimal_control/OCpy/test/oc/genetic/nanop/"
folder = "/Users/castd/Desktop/OCpy_hybrid_implementation/"

namefile = "input.dat"




OC_system=ini.SystemManager()
OC_system.init_system(folder, namefile)

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
