import SystemManager as ini
import time

folder = "/home/mana/Desktop/OC/cianidina/"
#folder = "/home/mana/programmi/python/optimal_control/OCpy/test/2order_prop/vac/"
namefile = "input_genetic.dat"




OC_system=ini.SystemManager()
OC_system.init_system(folder, namefile)



#start=time.time()
OC_system.oc.iterate()
#end=time.time()
#print("serial: " + str(end-start))

