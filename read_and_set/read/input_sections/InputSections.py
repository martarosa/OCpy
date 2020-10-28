import numpy as np

from read_and_set.read.input_sections.ABCSection import ABCSection


class SectionSystem(ABCSection):
    def __init__(self):
        super().__init__()
        self.section = 'SYSTEM'
        self.section_default_dictionary = {'folder'      : '',
                                           'name'        : 'output',
                                           'nstep'       : 'missing',
                                           'dt'          : 'missing',
                                           'oc_algorithm': 'none',
                                           'propagator'  : 'missing',
                                           'ibm_external_opt' : 'none'}
        self.section_dictionary = {}
        self.allowed_val = [['oc_algorithm', ['none', 'rabitzi', 'rabitzii', 'genetic', 'nelder-mead', 'bfgs', 'cg']],
                            ['propagator',   ['eulero_1order', 'eulero_2order', 'rabitz', 'quantum_trotter_suzuki']]]
        self.case_unsensitive_keys = ['oc_algorithm', 'propagator']


    def init_default_folder(self, folder):
        self.section_default_dictionary['folder'] = folder




class SectionField(ABCSection):
    def __init__(self):
        super().__init__()
        self.section = 'FIELD'
        self.section_default_dictionary = {'field_type': 'const',
                                           'fi'        : '0.01 0.01 0.01',
                                           'omega'     : '0 0 0',
                                           'sigma'     : '0',
                                           't0'        : '0',
                                           'name_field_file': 'none'}
        self.section_dictionary = {}
        self.allowed_val = [['field_type', ['const', 'pip', 'sin', 'gau', 'sum', 'genetic', 'read']]]
        self.case_unsensitive_keys = ['field_type']


    def convert_string_coefficients(self, key_string_coeff):
        self.section_dictionary[key_string_coeff] = np.asarray(self.section_dictionary[key_string_coeff].split()).astype(float)
        rows = int(self.section_dictionary[key_string_coeff].shape[0] / 3)
        self.section_dictionary[key_string_coeff].reshape(rows, 3)



class SectionWaveFunction(ABCSection):
    def __init__(self):
        super().__init__()
        self.section = 'WAVEFUNCTION'
        self.section_default_dictionary = {'name_ci' : 'ci_ini.inp',
                                           'name_ei' : 'ci_energy.inp',
                                           'name_mut': 'ci_mut.inp'}
        self.section_dictionary = {}
        self.allowed_val = []
        self.case_unsensitive_keys = []



class SectionMedium(ABCSection):
    def __init__(self):
        super().__init__()
        self.section = 'MEDIUM'
        self.section_default_dictionary ={ 'medium'          : 'vac',
                                           'name_vij'        : 'ci_pot.inp',
                                           'name_file_cavity': 'cavity.inp',
                                           'name_mat_SD'     :'mat_SD.inp',
                                           'name_q_reaction_field_dyn': 'np_bem.mdy',
                                           'name_q_local_field_dyn': 'np_bem.mld'}

        self.section_dictionary = {}
        self.allowed_val = [['medium', ['vac', 'sol', 'nanop']]]
        self.case_unsensitive_keys = ['medium', 'read_qijn']



class SectionSave(ABCSection):
    def __init__(self):
        super().__init__()
        self.section = "SAVE"
        self.section_default_dictionary = {'restart_step': '10', 'append': 'false'}
        self.section_dictionary = {}
        self.allowed_val = [['append', ['true', 'false']]]
        self.case_unsensitive_keys = ['append']


class SectionOptimalControl(ABCSection):
    def __init__(self):
        super().__init__()
        self.section = "OPTIMALC"
        self.section_default_dictionary = { 'restart': 'false',
                                            'alpha'  : 'const',
                                            'alpha_file': 'none',
                                            'target_state': '1',
                                            'n_iterations': '0',
                                            'convergence_thr': '99999',
                                            'conf_file': 'none'}
        self.section_dictionary = {}
        self.allowed_val = [['alpha', ['const', 'sin', 'quin']],
                            ['restart', ['true', 'false']]]
        self.case_unsensitive_keys = ['restart', 'alpha']



