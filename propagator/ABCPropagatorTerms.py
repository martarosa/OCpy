from abc import ABCMeta, abstractmethod

class ABCPropagatorTerms(metaclass=ABCMeta):
    def __init__(self):
        self.dict_terms = {}

    @abstractmethod
    def init(self, *args):
        pass
