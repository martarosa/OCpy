import SystemManager as ini
from deap import base
from deap import creator
from deap import tools

import random

import argparse
import time as time
#from DeapWrapper import DeapCreator


#python3 run_oc.py -f folder/namefile takes as input the string, cut at the last "/" and split folder and namefile
#if only the namefile is provided default folder is ./
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

folder = "/home/mana/programmi/python/optimal_control/OCpy/test/genetic/"
namefile = "input.dat"


#creator.create("J", base.Fitness, weights=(1.0,))
#creator.create("Chromosome", list, fitness=creator.J)
#toolbox = base.Toolbox()
#toolbox.register("create_random_ampl", random.uniform, -0.001, 0.001)
#toolbox.register("single_chromosome", tools.initIterate, creator.Chromosome, toolbox.create_random_ampl, n=1)








OC_system=ini.SystemManager()
OC_system.init_system(folder, namefile)

#deap_creator = DeapCreator()
#deap_creator.create()

#psi=system.wavef.propagate_n_step_2order(system.oc.dt,
#                                         system.oc.nstep,
#                                         system.starting_field.get_field(),
#                                         system.env.get_qijn(),
#                                         system.env.get_Vij())

#chi=system.wavef.propagate_one_step_bwd(system.oc.dt,system.)

#m=af.population_from_wf_matrix(psi)
#np.savetxt(folder+"propagate_pop.dat", np.real(m[:,0]),fmt='%0.4e')
#np.savetxt(folder+"pop_propagation2Eulero.dat", np.real(m))



#start=time.time()
#OC_system.oc.iterate()
#end=time.time()
#print(end-start)
#import auxiliary_functions as af

#doublesum = af.double_summation(system.wavef.ci,system.wavef.ci,system.env.get_qijn())

#single_sum_tessere = af.single_summation_tessere(doublesum, system.env.get_Vij())

#single_sum = af.single_summation(system.wavef.ci, single_sum_tessere)
