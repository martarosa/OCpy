from abc import ABCMeta, abstractmethod


class PCMParameters():
    def __init__(self):
        self.env = None
        self.cavity = None
        self.muLF = 0


class ABCPCM(metaclass=ABCMeta):
    def __init__(self):
        self.par = PCMParameters()

        self.q_t = None
        self.qijn = None
        self.qijn_lf = None
        self.q00n = None

    @abstractmethod
    def init_pcm(self, PCM_input, mol, field_t):
        pass

    @abstractmethod
    def propagate(self, i, mol, field_t):
        pass

    @abstractmethod
    def get_q_t(self):
        pass









