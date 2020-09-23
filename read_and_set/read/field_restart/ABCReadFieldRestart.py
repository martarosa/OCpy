from read_and_set.input.FieldInput import FieldInput
from abc import ABCMeta, abstractmethod

class ABCReadFieldRestart(metaclass=ABCMeta):
    def __init__(self):
        self.input_parameters = FieldInput()

    @abstractmethod
    def read_file(self, folder, namefile):
        pass



