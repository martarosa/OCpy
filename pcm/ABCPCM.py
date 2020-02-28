from abc import ABCMeta, abstractmethod

from parameters.PCMParameters import PCMParameters


class ABCPCM(metaclass=ABCMeta):
    def __init__(self):
        self.par = PCMParameters()

    @abstractmethod
    def init_pcm(self, PCM_input, mol, field_dt_vector):
        pass

    @abstractmethod
    def propagate(self, i, mol, field_dt_vector):
        pass

    @abstractmethod
    def get_q_t(self):
        pass




