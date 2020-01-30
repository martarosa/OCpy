import numpy as np

from field.FieldInput import FieldInput
from field.ReadFieldRestart import ReadFieldRestart
import configparser


class ReadFieldRestartGenetic(ReadFieldRestart):
    def __init__(self):
        super().__init__()
        self.field_par = FieldInput()


    def read_file(self, folder, namefile):
        load = np.loadtxt(folder + namefile) # w0, w1, ...wn \n a0 a1 ...an
        n = load.shape[0]/2
        self.field_par.omega = load[:n]
        self.field_par.fi = load[n:]
        self.field_par. field_type = 'restart_genetic'
        self.field_par.sigma = 0
        self.field_par.t0 = 0
        self.field_par.namefile = 'false'
        self.field_par.omega_max = 0
