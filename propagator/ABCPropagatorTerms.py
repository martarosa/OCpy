from abc import ABCMeta, abstractmethod

class ABCPropagatorTerms():
    def __init__(self):
        self.dict_terms = {}

    @abstractmethod
    def init(self, *args):
        pass
