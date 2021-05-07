import argparse
import SystemManager as ini
import time
from read_and_set.read.ReadOutputQuantumCalc import ReadOutputQuantumCalc


folder = "/home/mana/programmi/python/optimal_control/OCpy/test/oc/genetic/vac/2levels/"
namefile = "input.dat"


variabile = 5



OC_system=ini.SystemManager()
OC_system.init_system(folder, namefile)

#print(OC_system.mol.par.muT[2,7])
#print(OC_system.mol.par.muT[7,1])




#start=time.time()
OC_system.oc.iterate()
#end=time.time()
#print("serial: " + str(end-start))


#read = ReadOutputGaussian()
#read.convertVgamess_to_VWaveT(folder + "ci_pot.inp", folder + "ci_pot_transf.inp", 11)
