import numpy as np

from read_and_set.input.FieldInput import FieldInput
from read_and_set.read.field_restart.ABCReadFieldRestart import ABCReadFieldRestart


class ReadGeneticFieldRestart(ABCReadFieldRestart):
    def __init__(self):
        super().__init__()
        self.input_parameters = FieldInput()


    def read_file(self, folder, namefile):
        load = np.loadtxt(folder + namefile) # w0, w1, ...wn \n a0 a1 ...an
        n = int((load.shape[0])/2)
        self.input_parameters.fi = load[n:]
        self.input_parameters.omega = load[:n]
        self.input_parameters.field_type = 'restart_genetic'
        self.input_parameters.sigma = 0
        self.input_parameters.t0 = 0
        self.input_parameters.namefile = 'false'
        self.input_parameters.omega_max = 0
