import SystemManager as ini

folder = "/home/mana/programmi/python/optimal_control/OCpy/test/oc/rabitz/fortran_pcm_prop/"
namefile = "input.dat"




OC_system=ini.SystemManager()
OC_system.init_system(folder, namefile)



#start=time.time()
OC_system.oc.iterate()
#end=time.time()
#print(end-start)

