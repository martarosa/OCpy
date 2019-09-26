import SystemManager as ini
import argparse
import time as time
from deap import base, creator
import numpy as np
import genetic as oc
folder = "/home/mana/Desktop/python/optimal_control/OCpy/nanop/test/oc/fullpy/"
namefile = "input.dat"


OC_system=ini.SystemManager()
OC_system.init_system(folder, namefile)



gen_mio = oc.GeneticParameters()
gen_mio.set_parameters(0.01, 48, OC_system.oc.prop_psi.mol)







creator.create("Fitness", base.Fitness, weights=(-1.0,))