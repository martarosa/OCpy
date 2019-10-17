import numpy as np

from field.FieldParameters import FieldParameters
from field.ReadFieldRestart import ReadFieldRestart


class ReadFieldRestartRabitz(ReadFieldRestart):
    def __init__(self):
        super().__init__()
        self.field_par = FieldParameters()


    def read_file(self, folder, namefile):
        self.field_par.field = np.loadtxt(folder + namefile, usecols=(0, 1, 2))
        self.field_par. field_type = 'optimizedRabitz'
        self.field_par.fi = np.array([0, 0, 0])
        self.field_par.fi_cos = np.array([0, 0, 0])
        self.field_par.omega = np.array([0, 0, 0])
        self.field_par.sigma = 0
        self.field_par.t0 = 0
        self.field_par.namefile = 'false'
        self.field_par.omega_max = 0