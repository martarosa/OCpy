import numpy as np

from read_and_set.input.FieldInput import FieldInput
from read_and_set.read.ABCReadFieldRestart import ABCReadFieldRestart


class ReadFieldRestartRabitz(ABCReadFieldRestart):
    def __init__(self):
        super().__init__()
        self.input_parameters = FieldInput()


    def read_file(self, folder, namefile):
        self.input_parameters.field = np.loadtxt(folder + namefile, usecols=(1, 2, 3))
        self.input_parameters.field_time_axis = np.loadtxt(folder + namefile, usecols=(0))
        self.input_parameters.field_type = 'restart_rabitz'
        self.input_parameters.fi = None
        self.input_parameters.omega = None
        self.input_parameters.sigma = None
        self.input_parameters.t0 = None

        self.input_parameters.omega_sys = None
