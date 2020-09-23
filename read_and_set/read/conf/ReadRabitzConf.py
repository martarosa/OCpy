import configparser
import os.path
from copy import deepcopy

import read_and_set.read.RabitzConfSections

from read_and_set.read import auxiliary_functions as af
from read_and_set.read.input_sections.ABCReadInputFile import ABCReadInputFile


class ReadRabitzConf(ABCReadInputFile):
    def __init__(self):
        super().__init__()
        self.n_sections = None
        self.sections = ["PARAMETERS"]
        self.parameters = read_and_set.read.input_sections.RabitzConfSections.SectionRabitzParameters()
        self.conf_str = None


    def read_file(self, folder, namefile):
        if(not os.path.isfile(folder+namefile)):
            af.exit_error("ERROR: input file not found")
        user_input = configparser.ConfigParser()
        self.n_sections = len(self.sections)
        user_input.read(folder + namefile)
        self.set_nml_sections(user_input)
        self.check_input_sections_names(user_input)



    def set_nml_sections(self, user_input):
        self.parameters.init(user_input)
        self.create_print_string()


    def create_print_string(self):
        self.conf_str = deepcopy(self.parameters.section_dictionary)
        self.conf_str = str(self.conf_str)


