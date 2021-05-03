import configparser
import os.path
from read_and_set.read import auxiliary_functions as af
from read_and_set.read import NameListSections as sec
from read_and_set.read.ABCReadNamelist import ABCReadNamelist


class ReadNamelistOC(ABCReadNamelist):
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
        pass

    def check_wavef_nml_consistency(self):
        pass

    def check_field_nml_consistency(self):
        #if calculation is restarted only "name_field_file" can be present in FIELD nml
        if(self.oc.check_namelist_key_exist_and_value('restart', 'true')):
        #number of key is max equal to 1
            if (self.field.check_namelist_key_exist("fi") or 
               self.field.check_namelist_key_exist("sigma") or
               self.field.check_namelist_key_exist("t0") or 
               self.field.check_namelist_key_exist("omega") or
               self.field.check_namelist_key_exist("field_type")): 
                af.exit_error("ERROR. Restarted calculation. field is read from file. "
                                 "field parameters in namelist \"FIELD\" are not used ")
            self.check_field_restart()
        #if calculaion is with genetic oc algorithm only genetic field is allowed
        elif(self.sys.check_namelist_key_exist_and_value('oc_algorithm', 'genetic')):
                #if oc optimizator genetic only options are genetic field or restart, all other keywords in FIELD must be empty
                if (not self.field.check_namelist_key_exist_and_value('field_type', 'genetic')):
                    af.exit_error(
                        "Error. Genetic algorithm for optimal control optimization. "
                        "\"field_type\" in \"FIELD\" namelist must be \"genetic\"")
        self.check_field_shape_parameters_warning()

    def check_field_restart(self):
        #if there is a name for the restaring field check if it exist
        if self.field.check_namelist_key_exist('name_field_file'):
            if os.path.isfile(self.sys.section_dictionary['folder'] + self.field.section_dictionary['name_field_file']):
                pass
            #if it does not exist  checks if default name exist and make WARNiNG
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
        if (self.sys.check_namelist_key_exist_and_list_value('oc_algorithm', ['eulero_1order_prop', 'eulero_2order_prop'])):
            if len(self.oc.section_dictionary.keys()) != 0:
                af.exit_error("ERROR. \"oc_algorithm\" in namelist \"SYSTEM\" is a propagation without optimal control. "
                         "OPTIMALC namelist should be empty")

    def check_save_nml_consistency(self):
        pass



# set dependent default values
    def set_oc_dependent_default(self):
        #set default name of configuration file for oc iterator depending on algorithm
        if not self.oc.check_namelist_key_exist('iterator_config_file'):
            if self.sys.check_namelist_key_exist_and_value('oc_algorithm', 'genetic'):

                self.oc.add_key_and_value_to_namelist('iterator_config_file', 'genetic.conf')


    def set_save_dependent_default(self):
        # set default value of append = true if restart asked
        if (self.oc.check_namelist_key_exist_and_value('restart', 'true')):
            if not(self.save.check_namelist_key_exist('append')):
                self.save.add_key_and_value_to_namelist('append', 'true')



