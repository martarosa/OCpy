from abc import ABCMeta, abstractmethod


class PCMParameters():
    def __init__(self):
        self.env = None
        self.cavity = None
        self.muLF = 0


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









