from abc import ABCMeta
import configparser
import os.path
import numpy as np

import dictionaries.OCDictionaries

from read_and_set.read.input_sections import InputSections as sec
from read_and_set.read import auxiliary_functions as af
from read_and_set.read.input_sections.ABCReadInputFile import ABCReadInputFile

# Read input.dat were there are "SYSTEM" , "WAVEFUNCTION", "FIELD", "MEDIUM", "OPTIMALC", "SAVE" namelist
# even empty they must all be present
# For each of them
# 1-checks if required information are missing
# 2-lowers non case unsensitive keys
# 3-check allowed values for some keys
#Then performs other checks to give warning messages
#This file is commented but not polished neither organized

class ReadInput(ABCReadInputFile):
    def __init__(self):
        super().__init__()
        self.n_sections = None
        self.default_dict = None
        self.sections = ["SYSTEM" , "WAVEFUNCTION", "FIELD", "MEDIUM", "OPTIMALC", "SAVE"]
        self.sys = sec.SectionSystem()
        self.wf = sec.SectionWaveFunction()
        self.field = sec.SectionField()
        self.medium = sec.SectionMedium()
        self.save = sec.SectionSave()
        self.oc = sec.SectionOptimalControl()

    def read_file(self, folder, namefile):
        if(not os.path.isfile(folder+namefile)):
            af.exit_error("ERROR: input file not found")
        user_input = configparser.ConfigParser()
        self.n_sections = len(self.sections)
        user_input.read(folder + namefile)
        self.sys.init_default_folder(folder)
        self.set_nml_sections(user_input)
        self.check_input_sections_names(user_input)


    def set_nml_sections(self, user_input):
        self.sys.init(user_input)
        self.sys.check()
        self.field.init(user_input)
        self.field.check()
        self.wf.init(user_input)
        self.wf.check()
        self.save.init(user_input)
        self.save.check()
        self.medium.init(user_input)
        self.medium.check()
        self.oc.init(user_input)
        self.oc.check()

        self.check_system_nml_consistency()
        self.sys.fill_empty_with_default()

        self.check_field_nml_consistency()
        self.field.fill_empty_with_default()
        self.field.convert_string_coefficients('fi')
        self.field.convert_string_coefficients('omega')

        self.check_wavef_nml_consistency()
        self.wf.fill_empty_with_default()

        self.check_medium_nml_consistency()
        self.medium.fill_empty_with_default()

        self.check_optimalc_nml_consistency()
        self.set_oc_dependent_default()
        self.oc.fill_empty_with_default()

        self.check_save_nml_consistency()
        self.set_save_dependent_default()
        self.save.fill_empty_with_default()


        self.sys.check_keys()
        self.wf.check_keys()
        self.medium.check_keys()
        self.field.check_keys()
        self.oc.check_keys()
        self.save.check_keys()





    def check_system_nml_consistency(self):
        check = OCvsPropagatorCombinedKeywordValues()
        if not check.check(self.sys.section_dictionary["oc_algorithm"], self.sys.section_dictionary["propagator"]):
            af.exit_error("ERROR. \"propagator\" and \"oc_algorithm\" keyword are combined in a way that is not allowed")


    def check_wavef_nml_consistency(self):
        pass


    def check_field_nml_consistency(self):
        #if calculation is restarted only "name_field_file" can be present in FIELD nml
        if(self.oc.check_namelist_key_exist_and_value('restart', 'true')):
        #number of key is max equal to 1
            if len(self.field.section_dictionary.keys()) > 1:
                af.exit_error("ERROR. Restarted calculation. field is read from file. "
                                 "Keys in namelist \"FIELD\" are not used apart \"name_field_file\"")
        #and that one must be "name_field_file"
            elif len(self.field.section_dictionary.keys()) == 1:
                if self.field.check_namelist_key_exist("name_field_file"):
                    pass
                else:
                    af.exit_error(
                            "Error. Restarted calculation. "
                            "Field is read from file. Keys in namelist \"FIELD\" are not used apart \"name_field_file\"")
            self.check_field_restart()
        #checking consistenci with oc_algorithm
        else:
            check = OCvsFieldCombinedKeywordValues()
            if not check.check(self.sys.section_dictionary["oc_algorithm"], self.field.section_dictionary["field_type"]):
                af.exit_error("ERROR. \"field_type\" and \"oc_algorithm\" keyword are combined in a way that is not allowed")

        self.check_field_shape_parameters_warning()


    def check_field_restart(self):
        #if there is a name for the restaring field check if it exist
        if self.field.check_namelist_key_exist('name_field_file'):
            if os.path.isfile(self.sys.section_dictionary['folder'] + self.field.section_dictionary['name_field_file']):
                pass
            #if it does not exist  checks if default name exist and make WARNING
            elif os.path.isfile(self.sys.section_dictionary['folder'] + self.sys.section_dictionary['name'] + '_field_bkp.dat'):
                print("WARNING: restart from \"" + self.field.section_dictionary['name_field_file'] + "\" asked but file not found. \n"
                      "Restarting from " + self.sys.section_dictionary['name'] + '_field_bkp.dat')
                self.field.section_dictionary['name_field_file'] = self.sys.section_dictionary['name'] + '_field_bkp.dat'
                self.oc.section_dictionary['restart'] = "only_bkp_found"
            # if a restart don't exist start from new and makes a WARNING
            else:
                self.oc.section_dictionary['restart'] = 'norestart_found'
                print("WARNING: restart from \"" + self.field.section_dictionary['name_field_file'] + "\" asked but file not found. \n"
                      "Starting calculation from default field")
                self.field.section_dictionary['name_field_file'] = 'false'
        #f there is no name directly checks if default restart exist
        elif os.path.isfile(self.sys.section_dictionary['folder'] + self.sys.section_dictionary['name'] + '_field_bkp.dat'):
                self.field.section_dictionary['name_field_file'] = self.sys.section_dictionary['name'] + '_field_bkp.dat'
        # if don't find a restart start from new field and make a WARNING
        else:
            self.oc.section_dictionary['restart'] = 'norestart_found'
            print("WARNING: restart asked but file not found. \n"
                      "Starting calculation from default field")
            self.field.section_dictionary['name_field_file'] = 'false'

    def check_field_shape_parameters_warning(self):
        # if default (which is constant) or constant
        if(not self.field.check_namelist_key_exist('field_type')
                or self.field.check_namelist_key_exist_and_value('field_type', 'const')):
            self.field.check_namelist_key_and_print('omega',
                                    "WARNING. in \"FIELD\" namelist \"omega\" keyword is not used")
            self.field.check_namelist_key_and_print('sigma',
                                    "WARNING. in \"FIELD\" namelist \"sigma\" keyword is not used")
            self.field.check_namelist_key_and_print('t0',
                                "WARNING. in \"FIELD\" namelist \"t0\" keyword is not used")
        elif(self.field.section_dictionary['field_type'] == 'gau'):
            self.field.check_namelist_key_and_print('omega',
                                    "WARNING. in \"FIELD\" namelist \"omega\" keyword is not used")
        elif(self.field.section_dictionary['field_type'] == 'sum'):
            self.field.check_namelist_key_and_print('sigma',
                                    "WARNING. in \"FIELD\" namelist \"sigma\" keyword is not used")
            self.field.check_namelist_key_and_print('t0',
                                "WARNING. in \"FIELD\" namelist \"t0\" keyword is not used")
        elif(self.field.section_dictionary['field_type'] == 'genetic'):
            self.field.check_namelist_key_and_print('fi',
                                    "WARNING. in \"FIELD\" namelist \"fi\" keyword is not used, "
                                    "frequencies are assigned randomly")
            self.field.check_namelist_key_and_print('omega',
                                    "WARNING. in \"FIELD\" namelist \"omega\" keyword is not used, "
                                    "frequencies are chosen automatically depending on the genetic algorithm "
                                    "chosen")
            self.field.check_namelist_key_and_print('sigma',
                                    "WARNING. in \"FIELD\" namelist \"sigma\" keyword is not used")
            self.field.check_namelist_key_and_print('t0',
                                "WARNING. in \"FIELD\" namelist \"t0\" keyword is not used")

    def check_medium_nml_consistency(self):
        #if vac no names for medium files are allowed
        if (not self.medium.check_namelist_key_exist('medium')
                or self.medium.check_namelist_key_exist_and_value('medium', 'vac')):
            if len(self.medium.section_dictionary.keys()) > 1:
                af.exit_error("ERROR: in \"MEDIUM\" key \"medium\" is \"vac\". all other key in namelist MEDIUM are not used")

    def check_optimalc_nml_consistency(self):
        #if propagation OPTIMALC namelist must be empty
        if (self.sys.check_namelist_key_exist_and_list_value('oc_algorithm', ['none'])):
            if len(self.oc.section_dictionary.keys()) != 0:
                af.exit_error("ERROR. You are performing a propagation without optimization."
                         "OPTIMALC namelist should be empty")

    def check_save_nml_consistency(self):
        pass



# set dependent default values
    def set_oc_dependent_default(self):
        #set default name of configuration file for oc iterator depending on algorithm
        if not self.oc.check_namelist_key_exist('conf_file'):
            self.oc.add_key_and_value_to_namelist('conf_file',
                                                  dictionaries.OCDictionaries.OCConfigFileDefaultNames[self.sys.section_dictionary["oc_algorithm"]])


    def set_save_dependent_default(self):
        # set default value of append = true if restart asked
        if (self.oc.check_namelist_key_exist_and_value('restart', 'true')):
            if not(self.save.check_namelist_key_exist('append')):
                self.save.add_key_and_value_to_namelist('append', 'true')






class ABCCombinedKeywordValues(metaclass=ABCMeta):
    def __init__(self):
        self.list1 = []
        self.list2 = []
        self.matrix_allowed_couples = None


    def check(self, value1, value2):
        index1 = self.list1.index(value1)
        index2 = self.list2.index(value2)
        return self.matrix_allowed_couples[index1,index2]



class OCvsPropagatorCombinedKeywordValues(ABCCombinedKeywordValues):
    def __init__(self):
        super().__init__()

        self.list1 = ['none', 'rabitzi', 'rabitzii', 'genetic', 'nelder-mead', 'bfgs', 'cg']
        self.list2 = ['eulero_1order', 'eulero_2order', 'rabitz', 'quantum_trotter_suzuki']

        self.matrix_allowed_couples = np.array([[True, True, False, True],
                                                [False, False, True, False],
                                                [False, False, True, False],
                                                [True, True, False, True],
                                                [True, True, False, True],
                                                [True, True, False, True],
                                                [True, True, False, True]
                                               ])


class OCvsFieldCombinedKeywordValues(ABCCombinedKeywordValues):
    def __init__(self):
        super().__init__()

        self.list1 = ['none', 'rabitzi', 'rabitzii', 'genetic', 'nelder-mead', 'bfgs', 'cg']
        self.list2 = ['const', 'pip', 'sin', 'gau', 'sum', 'genetic', 'read']

        self.matrix_allowed_couples = np.array([[True, True, True, False, True, True, True],
                                                [True, True, True, False, True, True, True],
                                                [True, True, True, False, True, True, True],
                                                [True, True, True, False, True, True, True],
                                                [True, True, True, False, True, True, True],
                                                [False, False, False, False, False, False, False],
                                                [True, True, True, False, True, True, True]
                                                ])