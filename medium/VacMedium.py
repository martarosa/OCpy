from medium.ABCMedium import ABCMedium
from parameters.MediumParameters import MediumParameters

#as in vacuo there is no term for the medium, this class in empty and does nothing.
#It is initlized as it is needed in the propagation procedure to return 0 for the medium term

class VacMedium(ABCMedium):
    def __init__(self):
        self.par = MediumParameters()


    def init_medium(self, PCM_input, mol, field_dt_vector):
        self.par.cavity = 0
        self.par.medium = 'vac'

    def reset_medium(self, *args):
        pass

    def propagate_charges(self, mol, field_dt_vector):
        return 0
    def propagate_charges_fortran(self, mol, field_dt_vector):
        self.propagate(mol, field_dt_vector)


    def get_q_t(self):
        return 0
