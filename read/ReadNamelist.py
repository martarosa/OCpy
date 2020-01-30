import sys
from abc import ABCMeta, abstractmethod


class ReadNamelist(metaclass=ABCMeta):
    def __init__(self):
        self.n_sections = None
        self.default_dict = None
        self.sections = []

    @abstractmethod
    def set_default_dict(self, *args):
        pass

    @abstractmethod
    def read_file(self, folder, namefile):
        pass

    @abstractmethod
    def set_nml_sections(self, user_input):
        pass


    @abstractmethod
    def check_input_consistency(self, user_input):
        pass

    def check_input_sections_names(self, user_input):
        if not all(elem in user_input.sections() for elem in self.sections):
            sys.exit("Error. Sections names are wrong in input file")


    def check_namelist_key(self, key, output_string):
        if key in locals() or key in globals():
            print(output_string)
