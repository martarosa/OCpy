from abc import ABCMeta, abstractmethod

from parameters.MediumParameters import MediumParameters


class ABCMedium(metaclass=ABCMeta):
    def __init__(self):
        self.par = MediumParameters()

    @abstractmethod
    def init_medium(self, medium_input, mol, field_object):
        pass

    @abstractmethod
    def reset_medium(self, *args):
        pass

    @abstractmethod
    def propagate(self, mol, field_dt_vector):
        pass

    @abstractmethod
    def get_q_t(self):
        pass




