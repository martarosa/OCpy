import configparser
import sys
import NameListSections as sec

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
        user_input.read_dict({'SYSTEM': {
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
                                  'name_Vij': 'ci_pot.inp',
                                  'name_Q_tdplas': 'np_bem.mdy',
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
                              })
        self.n_sections = len(user_input.sections())
        user_input.read(folder + namefile)
        self.check_input_sections(user_input)
        self.set_dictionaries(user_input)




    def set_dictionaries(self, user_input):
        self.sys.par = dict(user_input.items("SYSTEM"))
        self.sys.lowering_case()
        self.sys.check_parameters_allowed_values()
        self.field.par = dict(user_input.items("FIELD"))
        self.field.lowering_case()
        self.field.check_parameters_allowed_values()
        self.field.convert_string_coefficients('fi')
        self.field.convert_string_coefficients('fi_cos')
        self.field.convert_string_coefficients('omega')
        self.wf.par = dict(user_input.items("WAVEFUNCTION"))
        self.wf.lowering_case()
        self.wf.check_parameters_allowed_values()
        self.save.par = dict(user_input.items("SAVE"))
        self.save.lowering_case()
        self.save.check_parameters_allowed_values()
        self.env.par = dict(user_input.items("ENVIRON"))
        self.env.lowering_case()
        self.env.check_parameters_allowed_values()
        self.oc.par = dict(user_input.items("OPTIMALC"))
        self.oc.lowering_case()
        self.oc.check_parameters_allowed_values()


    def check_input_sections(self, user_input):
        if len(user_input.sections()) != self.n_sections:
            sys.exit("Error. Sections names are wrong in input file")











