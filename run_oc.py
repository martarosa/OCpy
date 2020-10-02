import argparse
import SystemManager as ini
import time
from read_and_set.read.ReadOutputGaussian import ReadOutputGaussian
#parser = argparse.ArgumentParser()
#parser.add_argument('-f', action='store', dest='input')
#inputline=parser.parse_args()
#if inputline.input.find('/') == -1:
#    folder = "./"
#    namefile = inputline.input
#else:
#    split = inputline.input.rsplit('/',1)
#    folder = split[0]+"/"
#    namefile = split[1]




#folder = "/home/mana/programmi/python/optimal_control/OCpy/test_rabitz/"
folder = "/home/mana/programmi/python/optimal_control/OCpy/test/test_check_waveT_comparison/sol/"
namefile = "input.dat"




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
