import numpy as np
from read_and_set.read.ABCNameListSection import ABCNamelistSection


class SectionSystem(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = 'SYSTEM'
        self.section_default_dictionary = {'folder': '',
                                           'name': 'output',
                                           'dt': '0.01',
                                           'oc_algorithm': 'rabitzi'}
        self.section_dictionary = {}
        self.allowed_val = [['oc_algorithm', ['eulero_1order_prop', 'eulero_2order_prop', 'rabitzi', 'rabitzii', 'genetic']]]
        self.case_unsensitive_keys = ['oc_algorithm']


    def init_default_folder(self, folder):
        self.section_default_dictionary['folder'] = folder



class SectionField(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = 'FIELD'
        self.section_default_dictionary = {'field_type': 'const',
                                           'fi': '0.01 0.01 0.01',
                                           'omega': '0 0 0',
                                           'sigma': '0',
                                           't0': '0',
                                           'name_field_file': 'false',
                                           'nstep': '10000',
                                           'additional_steps': '0'}
        self.section_dictionary = {}
        self.allowed_val = [['field_type', ['const', 'pip', 'sin', 'gau', 'sum', 'sum_pip', 'genetic', 'test', 'read', 'read_genetic']]]
        self.case_unsensitive_keys = ['field_type']


    def convert_string_coefficients(self, key_string_coeff):
        self.section_dictionary[key_string_coeff] = np.asarray(self.section_dictionary[key_string_coeff].split()).astype(float)
        rows = int(self.section_dictionary[key_string_coeff].shape[0] / 3)
        self.section_dictionary[key_string_coeff].reshape(rows, 3)


class SectionWaveFunction(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = 'WAVEFUNCTION'
        self.section_default_dictionary = {'name_ci': 'ci_ini.inp',
                                           'name_ei': 'ci_energy.inp',
                                           'name_mut': 'ci_mut.inp'}
        self.section_dictionary = {}
        self.allowed_val = []
        self.case_unsensitive_keys = []



class SectionMedium(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = 'MEDIUM'
        self.section_default_dictionary ={ 'medium': 'vac',
                                           'polarization_charges': 'non-eq',
                                           'name_vij': 'ci_pot.inp',
                                           'name_file_cavity': 'cavity.inp',
                                           'name_mat_SD':'mat_SD.inp',
                                           'name_q_reaction_field_dyn': 'np_bem.mdy',
                                           'name_q_local_field_dyn': 'np_bem.mld'}

        self.section_dictionary = {}
        self.allowed_val = [['medium', ['vac', 'sol', 'nanop']]]
        self.case_unsensitive_keys = ['medium', 'read_qijn']



class SectionSave(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = "SAVE"
        self.section_default_dictionary = {'restart_step': '10', 'append': 'false'}
        self.section_dictionary = {}
        self.allowed_val = [['append', ['true', 'false']]]
        self.case_unsensitive_keys = ['append']


class SectionOptimalControl(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = "OPTIMALC"
        self.section_default_dictionary = { 'restart': 'false',
                                            'alpha': 'const',
                                            'alpha0': '1',
                                            'target_state': '1',
                                            'n_iterations': '0',
                                            'convergence_thr': '99999',
                                            'delta_ts': '0',
                                            'Ns': '0',
                                            'iterator_config_file': 'false'}
        self.section_dictionary = {}
        self.allowed_val = [['alpha', ['const', 'sin', 'quin']],
                            ['restart', ['true', 'false']]]
        self.case_unsensitive_keys = ['restart', 'alpha']





class SectionGenetic(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = "GENETIC"
        self.section_default_dictionary = {'n_chromosomes': '120',
                                           'n_selected_chr': '20',
                                           'genetic_algorithm' : 'sequential',
                                           'amplitude_lim': '0.05'}
        self.section_dictionary = {}
        self.allowed_val = [['genetic_algorithm', ['sequential', 'mixed']]]
        self.case_unsensitive_keys = ['genetic_algorithm']

class SectionMate(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = "MATE"
        self.section_default_dictionary = {'mate' : 'cxUniform',
                                           'mate_probability': '1'}
        self.section_dictionary = {}
        self.allowed_val = [['mate', ['cxUniform']]]
        self.case_unsensitive_keys = []


class SectionMutate(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = "MUTATE"
        self.section_default_dictionary = {'mutate':'mutGaussian',
                                           'mutate_probability': '0.2',
                                           'mutate_mu':'0',
                                           'mutate_starting_sigma': '0.01'}
        self.section_dictionary = {}
        self.allowed_val = [['mutate',['mutGaussian']]]
        self.case_unsensitive_keys = []


class SectionSelect(ABCNamelistSection):
    def __init__(self):
        super().__init__()
        self.section = "SELECT"
        self.section_default_dictionary = {'select':'selBest'}
        self.section_dictionary = {}
        self.allowed_val = [['select',['selBest']]]
        self.case_unsensitive_keys = []



