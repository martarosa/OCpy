import configparser
import sys
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


    def set_default_dict(self, folder):
        self.default_dict = {'SYSTEM': {
                                    'folder': folder,
                                    'name': 'output',
                                    'nstep': '10000',
                                    'dt': '0.01',
                                    'oc_algorithm': 'rabitzi'},
                              'WAVEFUNCTION': {
                                  'name_ci': 'ci_ini.inp',
                                  'name_ei': 'ci_energy.inp',
                                  'name_mut': 'ci_mut.inp'},
                              'FIELD': {
                                  'field_type': 'const',
                                  'fi': '0.01 0.01 0.01',
                                  'omega': '0 0 0',
                                  'sigma': '0',
                                  't0': '0',
                                  'name_field_file': 'false'},
                              'ENVIRON': {
                                  'env': 'vac',
                                  'name_vij': 'ci_pot.inp',
                                  'name_q_tdplas': 'np_bem.mdy',
                                  'read_qijn': 'false',
                                  'name_file_qijn': 'qijn.dat',
                                  'name_file_cavity': 'cavity.inp',
                                  'name_q_local_field': 'np_bem.mld'},
                              'OPTIMALC': {
                                  'restart': 'false',
                                  'alpha': 'const',
                                  'alpha0': '1',
                                  'target_state': '1',
                                  'n_iterations': '0',
                                  'convergence_thr': '99999',
                                  'delta_ts': '0',
                                  'Ns': '0',
                                  'iterator_config_file' : 'None'},
                              'SAVE': {'restart_step': '10'}
                              }


    def read_file(self, folder, namefile):
        user_input = configparser.ConfigParser()
        #first read without default values for checks on consistency
        user_input.read(folder + namefile)
        self.check_input_sections_names(user_input)
        self.check_input_consistency(user_input)
        #after checks on consistency we read again with default
        user_input.read_dict(self.default_dict)
        self.n_sections = len(user_input.sections())
        user_input.read(folder + namefile)
        self.set_nml_sections(user_input)



    def set_nml_sections(self, user_input):
        self.sys.init(user_input, self.default_dict, "SYSTEM")
        self.sys.check()
        self.field.init(user_input, self.default_dict, "FIELD")
        self.field.check()
        self.field.convert_string_coefficients('fi')
        self.field.convert_string_coefficients('omega')
        self.wf.init(user_input, self.default_dict, "WAVEFUNCTION")
        self.wf.check()
        self.save.init(user_input, self.default_dict, "SAVE")
        self.save.check()
        self.env.init(user_input, self.default_dict, "ENVIRON")
        self.env.check()
        self.oc.init(user_input, self.default_dict, "OPTIMALC")
        self.oc.check()



    def check_input_consistency(self, user_input):
        self.check_system_nml_consistency(user_input)
        self.check_wavef_nml_consistency(user_input)
        self.check_field_nml_consistency(user_input)
        self.check_environ_nml_consistency(user_input)
        self.check_optimalc_nml_consistency(user_input)
        self.check_save_nml_consistency(user_input)


    def check_system_nml_consistency(self, user_input):
        pass

    def check_wavef_nml_consistency(self, user_input):
        pass





    def check_field_nml_consistency(self, user_input):
        #if restart only name_field_file can be present in FIELD nml
        if(user_input['OPTIMALC']['restart'].lower() == 'true'):
            if len(user_input['FIELD'].keys()) > 1:
                sys.exit("ERROR. Restarted calculation. field is read from file. "
                         "Keys in namelist \"FIELD\" are not used apart \"name_field_file\"")
            elif len(user_input['FIELD'].keys()) == 1:
                try: user_input['FIELD']["name_field_file"]
                except NameError:
                    sys.exit(
                        "Error. Restarted calculation. "
                        "Field is read from file. Keys in namelist \"FIELD\" are not used apart \"name_field_file\"")
        else:
            if (user_input['SYSTEM']['oc_algorithm'].lower() == 'genetic'):
                #if oc optimizator genetic only options are genetic field or restart, all other keywords in FIELD must be empty
                if (user_input['FIELD']['field_type'].lower() != 'genetic'):
                    sys.exit(
                        "Error. Genetic algorithm for optimal control optimization. "
                        "\"field_type\" in \"FIELD\" namelist must be \"genetic\"")
                if len(user_input['FIELD'].keys()) > 1:
                    sys.exit("Error. Genetic algorithm for optimal control optimization. "
                             "\"FIELD\" namelist should be empty apart for  \"field_type\"")
            else:
                self.check_field_shape_parameters_warning(user_input)

    def check_field_shape_parameters_warning(self, user_input):
        if(user_input['FIELD']['field_type'].lower() == 'const'):
            self.check_namelist_key(user_input['FIELD']['omega'],
                                    "WARNING. in \"FIELD\" namelist \"omega\" keyword is not used")
            self.check_namelist_key(user_input['FIELD']['sigma'],
                                    "WARNING. in \"FIELD\" namelist \"sigma\" keyword is not used")
            self.check_namelist_key(user_input['FIELD']['t0'],
                                "WARNING. in \"FIELD\" namelist \"t0\" keyword is not used")
        elif(user_input['FIELD']['field_type'].lower() == 'gau'):
            self.check_namelist_key(user_input['FIELD']['omega'],
                                    "WARNING. in \"FIELD\" namelist \"omega\" keyword is not used")
        elif(user_input['FIELD']['field_type'].lower() == 'sum'):
            self.check_namelist_key(user_input['FIELD']['sigma'],
                                    "WARNING. in \"FIELD\" namelist \"sigma\" keyword is not used")
            self.check_namelist_key(user_input['FIELD']['t0'].lower(),
                                "WARNING. in \"FIELD\" namelist \"t0\" keyword is not used")
        if(user_input['FIELD']['field_type'].lower() == 'genetic'):
            self.check_namelist_key(user_input['FIELD']['fi'],
                                    "WARNING. in \"FIELD\" namelist \"fi\" keyword is not used, "
                                    "frequencies are assigned randomly")
            self.check_namelist_key(user_input['FIELD']['omega'],
                                    "WARNING. in \"FIELD\" namelist \"omega\" keyword is not used, "
                                    "frequencies are chosen automatically depending on the genetic algorithm "
                                    "chosen")
            self.check_namelist_key(user_input['FIELD']['sigma'],
                                    "WARNING. in \"FIELD\" namelist \"sigma\" keyword is not used")
            self.check_namelist_key(user_input['FIELD']['t0'],
                                "WARNING. in \"FIELD\" namelist \"t0\" keyword is not used")

    def check_environ_nml_consistency(self, user_input):
        #if vac no names for pcm files are allowed
        if user_input['ENVIRON']['env'].lower() == 'vac':
            if len(user_input['ENVIRON'].keys()) > 1:
                sys.exit("ERROR: in \"ENVIRON\" key \"env\" is \"vac\". all other key in namelist ENVIRON are not used")

    def check_optimalc_nml_consistency(self, user_input):
        #if propagation OPTIMALC namelist mus be empty
        if (user_input['SYSTEM']['oc_algorithm'].lower() == 'eulero_1order_prop' or
                user_input['SYSTEM']['oc_algorithm'].lower() == 'eulero_2order_prop'):
            if len(user_input['OPTIMALC'].keys()) != 0:
                sys.exit("ERROR. \"oc_algorithm\" in namelist \"SYSTEM\" is a propagation without optimal control. "
                         "OPTIMALC namelist should be empty")

    def check_save_nml_consistency(self, user_input):
        pass






class ReadNamelistOC_old():
    def __init__(self):
        self.n_sections = None
        self.sys = sec.SectionSystem()
        self.wf = sec.SectionWaveFunction()
        self.field = sec.SectionField()
        self.env = sec.SectionEnviron()
        self.save = sec.SectionSave()
        self.oc = sec.SectionOptimalControl()


    def read_file(self, folder, namefile):
        user_input = configparser.ConfigParser()
        default_dict = {'SYSTEM': {
                                    'folder': folder,
                                    'name': 'output',
                                    'nstep': '10000',
                                    'dt': '0.01',
                                    'oc_algorithm': 'rabitzi'},
                              'WAVEFUNCTION': {
                                  'name_ci': 'ci_ini.inp',
                                  'name_ei': 'ci_energy.inp',
                                  'name_mut': 'ci_mut.inp'},
                              'FIELD': {
                                  'field_type': 'const',
                                  'fi': '0.01 0.01 0.01',
                                  'omega': '0 0 0',
                                  'sigma': '0',
                                  't0': '0',
                                  'name_field_file': 'false'},
                              'ENVIRON': {
                                  'env': 'vac',
                                  'name_vij': 'ci_pot.inp',
                                  'name_q_tdplas': 'np_bem.mdy',
                                  'read_qijn': 'false',
                                  'name_file_qijn': 'qijn.dat',
                                  'name_file_cavity': 'cavity.inp',
                                  'name_q_local_field': 'np_bem.mld'},
                              'OPTIMALC': {
                                  'oc_algorithm': 'rabitzi',
                                  'restart': 'false',
                                  'alpha': 'const',
                                  'alpha0': '1',
                                  'target_state': '1',
                                  'n_iterations': '0',
                                  'convergence_thr': '99999',
                                  'delta_ts': '0',
                                  'Ns': '0',
                                  'iterator_config_file' : 'None'},
                              'SAVE': {'restart_step': '10'}
                              }
        #first read without default values for checks on consistency
        user_input.read(folder + namefile)
        self.check_input_sections(user_input)
        #after checks on consistency we read again with default
        user_input.read_dict(default_dict)
        self.n_sections = len(user_input.sections())
        user_input.read(folder + namefile)
        self.set_dictionaries(user_input, default_dict)



    def set_dictionaries(self, user_input, default_dict):
        self.sys.init(user_input, default_dict, "SYSTEM")
        self.sys.check()
        self.field.init(user_input, default_dict, "FIELD")
        self.field.check()
        self.field.convert_string_coefficients('fi')
        self.field.convert_string_coefficients('omega')
        self.wf.init(user_input, default_dict, "WAVEFUNCTION")
        self.wf.check()
        self.save.init(user_input, default_dict, "SAVE")
        self.save.check()
        self.env.init(user_input, default_dict, "ENVIRON")
        self.env.check()

        self.oc.section_dictionary = dict(user_input.items("OPTIMALC"))
        self.oc.check()


    def check_input_sections(self, user_input):
        # if len(user_input.sections()) != self.n_sections:
        sections = ["SYSTEM" , "WAVEFUNCTION", "FIELD", "ENVIRON", "OPTIMALC", "SAVE"]
        if not all(elem in user_input.sections() for elem in sections):
            sys.exit("Error. Sections names are wrong in input file")









