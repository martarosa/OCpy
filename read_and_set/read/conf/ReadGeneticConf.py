import configparser
import os.path
from copy import deepcopy

import read_and_set.read.GeneticConfSections

from read_and_set.read import auxiliary_functions as af
from read_and_set.read.input_sections.ABCReadInputFile import ABCReadInputFile


class ReadGeneticConf(ABCReadInputFile):
    def __init__(self):
        super().__init__()
        self.n_sections = None
        self.sections = ["GENETIC", "MATE", "MUTATE", "SELECT"]
        self.genetic = read_and_set.read.input_sections.GeneticConfSections.SectionGenetic()
        self.mate = read_and_set.read.input_sections.GeneticConfSections.SectionMate()
        self.mutate = read_and_set.read.input_sections.GeneticConfSections.SectionMutate()
        self.select = read_and_set.read.input_sections.GeneticConfSections.SectionSelect()
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
        self.genetic.init(user_input)
        self.genetic.check()
        self.mate.init(user_input)
        self.mate.check()
        self.mutate.init(user_input)
        self.mutate.check()
        self.select.init(user_input)
        self.select.check()

        self.check_genetic_nml_consistency()
        self.genetic.fill_empty_with_default()
        self.check_mate_nml_consistency()
        self.mate.fill_empty_with_default()
        self.check_mutate_nml_consistency()
        self.mutate.fill_empty_with_default()
        self.check_select_nml_consistency()
        self.select.fill_empty_with_default()

        self.genetic.check_keys()
        self.mate.check_keys()
        self.mutate.check_keys()
        self.select.check_keys()

        self.create_print_string()

    def create_print_string(self):
        self.conf_str = deepcopy(self.genetic.section_dictionary)
        self.conf_str.update(self.mate.section_dictionary)
        self.conf_str.update(self.mutate.section_dictionary)
        self.conf_str.update(self.select.section_dictionary)
        if self.conf_str['genetic_algorithm'] != 'mixed':
            del self.conf_str['eta_thr']
            del self.conf_str['q']
        self.conf_str = str(self.conf_str)



    def check_genetic_nml_consistency(self):
        pass

    def check_mate_nml_consistency(self):
        pass

    def check_mutate_nml_consistency(self):
        if(not self.genetic.check_namelist_key_exist("genetic_algorithm") or
                self.genetic.check_namelist_key_exist_and_value("genetic_algorithm", "sequential")):
            if(self.mutate.check_namelist_key_exist('eta_thr') or
                    self.mutate.check_namelist_key_exist('q')):
                print("WARNING: \"genetic_algorithm\" value equal \"sequential\""
                      "\"etha_thr\" and \"q\" keywords are not used")


    def check_select_nml_consistency(self):
        pass







