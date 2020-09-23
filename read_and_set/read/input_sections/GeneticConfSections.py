from read_and_set.read.input_sections.ABCSection import ABCSection


class SectionGenetic(ABCSection):
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


class SectionMate(ABCSection):
    def __init__(self):
        super().__init__()
        self.section = "MATE"
        self.section_default_dictionary = {'mate' : 'cxUniform',
                                           'mate_probability': '1'}
        self.section_dictionary = {}
        self.allowed_val = [['mate', ['cxUniform']]]
        self.case_unsensitive_keys = []


class SectionMutate(ABCSection):
    def __init__(self):
        super().__init__()
        self.section = "MUTATE"
        self.section_default_dictionary = {'mutate':'mutGaussian',
                                           'mutate_probability': '0.2',
                                           'n_mutate': '20',
                                           'mutate_mu':'0',
                                           'mutate_starting_sigma': '0.01',
                                           'eta_thr': '0.6',
                                           'q': '0.9'}
        self.section_dictionary = {}
        self.allowed_val = [['mutate',['mutGaussian']]]
        self.case_unsensitive_keys = []


class SectionSelect(ABCSection):
    def __init__(self):
        super().__init__()
        self.section = "SELECT"
        self.section_default_dictionary = {'select':'selBest'}
        self.section_dictionary = {}
        self.allowed_val = [['select',['selBest']]]
        self.case_unsensitive_keys = []