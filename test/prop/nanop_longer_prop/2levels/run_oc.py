import sys
sys.path.append('/scratch/martarosa/python/ocpy_tdplas_10_2020/ocpy/')
import argparse
import SystemManager as ini
import time
from read_and_set.read.ReadOutputGaussian import ReadOutputGaussian


folder = "./"
namefile = "input.dat"



OC_system=ini.SystemManager()
OC_system.init_system(folder, namefile)

OC_system.oc.iterate()
