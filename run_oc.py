import SystemManager as ini
import time
import sys
from multiprocessing import freeze_support

folder = "/Users/castd/Desktop/OCpy_restart_qiskit/"
namefile = "input.dat"


#import warnings

#warnings.simplefilter("ignore", DeprecationWarning)


variabile = 5

#if __name__ == "__main__":
    #freeze_support()
OC_system=ini.SystemManager()
OC_system.init_system(folder, namefile)
#start=time.time()
OC_system.oc.iterate()
#end=time.time()
#print("serial: " + str(end-start))

