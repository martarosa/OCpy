import configparser
import sys
import os.path
import read.auxiliary_functions as af
from read import NameListSections as sec
from read import ABCReadNamelist as readnamelist



class ReadNamelistOC(readnamelist.ABCReadNamelist):
    def __init__(self):
        super().__init__()
        self.n_sections = None
        self.default_dict = None
        self.sections = ["SYSTEM" , "WAVEFUNCTION", "FIELD", "ENVIRON", "OPTIMALC", "SAVE"]
        self.sys = sec.SectionSystem()
        self.wf = sec.SectionWaveFunction()
        self.field = sec.SectionField()
        self.env = sec.SectionEnviron()
        self.save = sec.SectionSave()
        self.oc = sec.SectionOptimalControl()


    def read_file(self, folder, namefile):
        if(not os.path.isfile(folder+namefile)):
            af.exit_error("ERROR: input file not found")
        user_input = configparser.ConfigParser()
        #first read without default values for checks on consistency
        self.n_sections = len(self.sections)
        user_input.read(folder + namefile)
        self.sys.init_default_folder(folder)
        self.set_nml_sections(user_input)


    def set_nml_sections(self, user_input):
        self.sys.init(user_input)
        self.sys.check()
        self.field.init(user_input)
        self.field.check()
        self.wf.init(user_input)
        self.wf.check()
        self.save.init(user_input)
        self.save.check()
        self.env.init(user_input)
        self.env.check()
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
        self.check_save_nml_consistency()
        self.save.fill_empty_with_default()
        self.check_environ_nml_consistency()
        self.env.fill_empty_with_default()
        self.check_optimalc_nml_consistency()
        self.oc.fill_empty_with_default()

        self.sys.check_keys()
        self.wf.check_keys()
        self.env.check_keys()
        self.field.check_keys()
        self.oc.check_keys()
        self.save.check_keys()

    def check_system_nml_consistency(self):
        pass

    def check_wavef_nml_consistency(self):
        pass

    def check_field_nml_consistency(self):
        #if restart only name_field_file can be present in FIELD nml
        #if(self.oc.section_dictionary('restart') in locals() or self.oc.section_dictionary('restart') in globals()):
        #    if( self.oc.section_dictionary('restart') == 'true'):
        if(self.oc.check_namelist_key_exist_and_value('restart', 'true')):
            if len(self.field.section_dictionary.keys()) > 1:
                af.exit_error("ERROR. Restarted calculation. field is read from file. "
                                 "Keys in namelist \"FIELD\" are not used apart \"name_field_file\"")

            elif len(self.field.section_dictionary.keys()) == 1:
                if self.field.check_namelist_key_exist("name_field_file"):
                    pass
                else:
                    af.exit_error(
                            "Error. Restarted calculation. "
                            "Field is read from file. Keys in namelist \"FIELD\" are not used apart \"name_field_file\"")
            self.check_field_resart()
        elif(self.sys.check_namelist_key_exist_and_value('oc_algorithm', 'genetic')):
                #if oc optimizator genetic only options are genetic field or restart, all other keywords in FIELD must be empty
                if (not self.field.check_namelist_key_exist_and_value('field_type', 'genetic')):
                    af.exit_error(
                        "Error. Genetic algorithm for optimal control optimization. "
                        "\"field_type\" in \"FIELD\" namelist must be \"genetic\"")
                if len(self.field.section_dictionary.keys()) > 1:
                    af.exit_error("Error. Genetic algorithm for optimal control optimization. "
                             "\"FIELD\" namelist should be empty apart for  \"field_type\"")
        self.check_field_shape_parameters_warning()


    def check_field_resart(self):
        if self.field.check_namelist_key_exist('name_field_file'):
            if os.path.isfile(self.sys.section_dictionary['folder'] + self.field.section_dictionary['name_field_file']):
                pass
            elif os.path.isfile(self.sys.section_dictionary['folder'] + self.sys.section_dictionary['name'] + '_field_bkp.dat'):
                print("WARNING: restart from \"" + self.field.section_dictionary['name_field_file'] + "\" asked but file not found. \n"
                      "Restarting from " + self.sys.section_dictionary['name'] + '_field_bkp.dat')
                self.field.section_dictionary['name_field_file'] = self.sys.section_dictionary['name'] + '_field_bkp.dat'
                self.oc.section_dictionary['restart'] = "only_bkp_found"
            else:
                self.oc.section_dictionary['restart'] = 'norestart_found'
                print("WARNING: restart from \"" + self.field.section_dictionary['name_field_file'] + "\" asked but file not found. \n"
                      "Starting calculation from default field")
                self.field.section_dictionary['name_field_file'] = 'false'
        elif os.path.isfile(self.sys.section_dictionary['folder'] + self.sys.section_dictionary['name'] + '_field_bkp.dat'):
                self.field.section_dictionary['name_field_file'] = self.sys.section_dictionary['name'] + '_field_bkp.dat'
        else:
            self.oc.section_dictionary['restart'] = 'norestart_found'
            print("WARNING: restart asked but file not found. \n"
                      "Starting calculation from default field")
            self.field.section_dictionary['name_field_file'] = 'false'




    def check_field_shape_parameters_warning(self):

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

    def check_environ_nml_consistency(self):
        #if vac no names for pcm files are allowed
        if (not self.env.check_namelist_key_exist('env')
                or self.env.check_namelist_key_exist_and_value('env', 'vac')):
            if len(self.env.section_dictionary.keys()) > 1:
                af.exit_error("ERROR: in \"ENVIRON\" key \"env\" is \"vac\". all other key in namelist ENVIRON are not used")

    def check_optimalc_nml_consistency(self):
        #if propagation OPTIMALC namelist mus be empty
        if (self.sys.check_namelist_key_exist_and_list_value('oc_algorithm', ['eulero_1order_prop', 'eulero_2order_prop'])):
            if len(self.oc.section_dictionary.keys()) != 0:
                af.exit_error("ERROR. \"oc_algorithm\" in namelist \"SYSTEM\" is a propagation without optimal control. "
                         "OPTIMALC namelist should be empty")



    def check_save_nml_consistency(self):
        pass




