import configparser
from read_and_set.read import NameListSections as sec
from read_and_set.read.ABCReadNamelist import ABCReadNamelist
import os.path
from read_and_set.read import auxiliary_functions as af
from copy import deepcopy

class ReadNamelistGenetic(ABCReadNamelist):
    def __init__(self):
        super().__init__()
        self.n_sections = None
        self.sections = ["GENETIC", "MATE", "MUTATE", "SELECT"]
        self.genetic = sec.SectionGenetic()
        self.mate = sec.SectionMate()
        self.mutate = sec.SectionMutate()
        self.select = sec.SectionSelect()
        self.string_file_config = None


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
        self.string_file_config = deepcopy(self.genetic.section_dictionary)
        self.string_file_config.update(self.mate.section_dictionary)
        self.string_file_config.update(self.mutate.section_dictionary)
        self.string_file_config.update(self.select.section_dictionary)
        if self.string_file_config['genetic_algorithm'] != 'mixed':
            del self.string_file_config['eta_thr']
            del self.string_file_config['q']
        self.string_file_config = str(self.string_file_config)



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







