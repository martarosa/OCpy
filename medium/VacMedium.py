from abc import ABCMeta, abstractmethod
from medium.ABCMedium import ABCMedium
from parameters.MediumParameters import MediumParameters


class VacMedium(ABCMedium):
    def __init__(self):
        self.par = MediumParameters()


    def init_medium(self, PCM_input, mol, field_dt_vector):
        self.par.cavity = 0
        self.par.medium = 'vac'


    def propagate(self, mol, field_dt_vector):
        return 0


    def get_q_t(self):
        return 0
