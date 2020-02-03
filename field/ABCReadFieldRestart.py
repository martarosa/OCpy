from field.FieldInput import FieldInput
from abc import ABCMeta, abstractmethod

class ABCReadFieldRestart(metaclass=ABCMeta):
    def __init__(self):
        self.field_par = FieldInput()

    @abstractmethod
    def read_file(self, folder, namefile):
        pass



