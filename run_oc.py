import SystemManager as ini
from read.ReadNamelistOC import ReadNamelistOC

folder = "/home/marta/programmi/python/optimal_control/OCpy/test/read_input/OC_input/nuovo/"

namefile = "input.dat"



user_input = ReadNamelistOC()
user_input.read_file(folder, namefile)



#OC_system=ini.SystemManager()
#OC_system.init_system(folder, namefile)



#start=time.time()
#OC_system.oc.iterate()
#end=time.time()
#print(end-start)

