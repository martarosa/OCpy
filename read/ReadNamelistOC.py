import configparser
import sys
from read import NameListSections as sec


class ReadNamelistOC():
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
                                    'propagation': 'eulero_2order'},
                              'WAVEFUNCTION': {
                                  'name_ci': 'ci_ini.inp',
                                  'name_ei': 'ci_energy.inp',
                                  'name_mut': 'ci_mut.inp'},
                              'FIELD': {
                                  'field_type': 'const',
                                  'fi': '0 0 0',
                                  'fi_cos':'0 0 0',
                                  'omega': '0 0 0',
                                  'sigma': '0',
                                  't0': '0',
                                  'name_field_file': 'false',
                                  'internal_check_genetic_field' : True  },
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
        user_input.read_dict(default_dict)
        self.n_sections = len(user_input.sections())
        user_input.read(folder + namefile)
        self.check_input_sections(user_input)
        self.set_dictionaries(user_input, default_dict)



    def set_dictionaries(self, user_input, default_dict):
        self.sys.init(user_input, default_dict, "SYSTEM")
        self.sys.check()
        self.field.init(user_input, default_dict, "FIELD")
        self.field.check()
        self.field.convert_string_coefficients('fi')
        self.field.convert_string_coefficients('fi_cos')
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








