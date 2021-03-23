import sys
sys.path.append('/home/mana/programmi/python/optimal_control/OCpy/')
import argparse
import SystemManager as ini
import time
from read_and_set.read.ReadOutputQuantumCalc import ReadOutputQuantumCalc


folder = "./"
namefile = "input.dat"



OC_system=ini.SystemManager()
OC_system.init_system(folder, namefile)

OC_system.oc.iterate()

