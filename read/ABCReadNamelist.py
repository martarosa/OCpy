import read.auxiliary_functions as af
from abc import ABCMeta, abstractmethod




class ABCReadNamelist(metaclass=ABCMeta):
    def __init__(self):
        self.n_sections = None
        self.sections = []

    @abstractmethod
    def read_file(self, folder, namefile):
        pass

    @abstractmethod
    def set_nml_sections(self, user_input):
        pass



    def check_input_sections_names(self, user_input):
        if not all(elem in user_input.sections() for elem in self.sections):
            af.exit_error("ERROR. Sections names are wrong in input file")



