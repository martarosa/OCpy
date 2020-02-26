from abc import ABCMeta, abstractmethod

from parameters.PCMParameters import PCMParameters


class VacPCM(metaclass=ABCMeta):
    def __init__(self):
        self.par = PCMParameters()


    def init_pcm(self, PCM_input, mol, field_dt_vector):
        self.par.muLF = 0
        self.par.cavity = 0
        self.par.env = 'vac'


    def propagate(self, i, mol, field_dt_vector):
        return 0


    def get_q_t(self):
        return 0