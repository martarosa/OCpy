import numpy as np

from field.FieldInput import FieldInput
from field.ABCReadFieldRestart import ABCReadFieldRestart


class ReadFieldRestartRabitz(ABCReadFieldRestart):
    def __init__(self):
        super().__init__()
        self.field_par = FieldInput()


    def read_file(self, folder, namefile):
        self.field_par.field = np.loadtxt(folder + namefile, usecols=(0, 1, 2))
        self.field_par. field_type = 'restart_rabitz'
        self.field_par.fi = np.array([0, 0, 0])
        self.field_par.omega = np.array([0, 0, 0])
        self.field_par.sigma = 0
        self.field_par.t0 = 0
        self.field_par.namefile = 'false'
        self.field_par.omega_max = 0
