import argparse
import SystemManager as ini
import time

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




#folder = "/home/mana/programmi/python/optimal_control/OCpy/test/tdplas-nanop/9-2020/"
folder = "/home/mana/programmi/python/optimal_control/OCpy/test/2order_prop/pcm/3sfere/"
namefile = "input.dat"




OC_system=ini.SystemManager()
OC_system.init_system(folder, namefile)



#start=time.time()
OC_system.oc.iterate()
#end=time.time()
#print("serial: " + str(end-start))

