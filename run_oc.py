import SystemManager as ini
from deap import base
from deap import creator
from deap import tools
from read import auxiliary_functions as af
import sys
sys.path.append('/home/mana/programmi/WaveT/TDPlas/src/')
from interface_tdplas import interface_tdplas as tdplas

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
import numpy as np
#folder = "/home/mana/programmi/python/optimal_control/OCpy/test/oc/fullpy/"
folder = "/home/mana/programmi/python/optimal_control/OCpy/test/nanoparticella/"
#folder ="/home/mana/programmi/python/optimal_control/OCpy/test/oc/"
namefile = "input.dat"



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

V_tdplas = af.double_summation(OC_system.mol.wf.ci_prev[0], np.conj(OC_system.mol.wf.ci_prev[0]), OC_system.mol.Vijn)

mut_tdplas = af.flip_3D_py2f(OC_system.mol.muT)


tdplas.set_tdplas(OC_system.oc.oc_iterator.class_attributes.dt,
                  OC_system.mol.wf.n_ci,
                  np.asfortranarray(OC_system.mol.wf.ci, dtype=np.complex64),
                  np.asfortranarray(OC_system.mol.en_ci, dtype=np.float64),
                  np.asfortranarray(mut_tdplas, dtype=np.float64), folder+"namelist_tdplas.inp")

#tdplas.init_medium(np.asfortranarray(V_tdplas, dtype=np.float64),
#                   np.asfortranarray(V_tdplas, dtype=np.float64))

#q = np.zeros(V_tdplas.size)
#tdplas.get_charge(np.asfortranarray(q, dtype=np.float64))




#start=time.time()
#OC_system.oc.iterate()
#end=time.time()
#print(end-start)
#import auxiliary_functions as af

#doublesum = af.double_summation(system.wavef.ci,system.wavef.ci,system.env.get_qijn())

#single_sum_tessere = af.single_summation_tessere(doublesum, system.env.get_Vij())

#single_sum = af.single_summation(system.wavef.ci, single_sum_tessere)
