from medium.ABCMedium import ABCMedium
from parameters.MediumParameters import MediumParameters


class VacMedium(ABCMedium):
    def __init__(self):
        self.par = MediumParameters()


    def init_medium(self, PCM_input, mol, field_dt_vector):
        self.par.cavity = 0
        self.par.medium = 'vac'

    def reset_medium(self, *args):
        pass

    def propagate(self, mol, field_dt_vector):
        return 0
    def propagate_fortran(self, mol, field_dt_vector):
        self.propagate(mol, field_dt_vector)


    def get_q_t(self):
        return 0