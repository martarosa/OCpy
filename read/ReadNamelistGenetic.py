import configparser
import sys
from read import NameListSections as sec


class ReadNamelistGenetic():
    def __init__(self):
        self.genetic = sec.SectionGenetic()


    def read_file(self, folder, namefile):
        user_input = configparser.ConfigParser()
        user_input.read_dict({'GENETIC': {
                                  'amplitude_min': '0.01',
                                  'amplitude_max': '0.01',
                                  'n_chromosomes': '1',
                                  'n_evolved_chr': 1}
                              })
        user_input.read(folder + namefile)
        self.check_input_sections(user_input)
        self.set_dictionaries(user_input)

    def set_dictionaries(self, user_input):
        self.genetic.par = dict(user_input.items("GENETIC"))
        self.genetic.lowering_case()
        self.genetic.check_parameters_allowed_values()


    def check_input_sections(self, user_input):
        if len(user_input.sections()) != 1:
            sys.exit("Error. Sections names are wrong in input file")