from read_and_set.set.ABCSetInput import ABCSetInput
from read_and_set.input.OCGeneticInput import OCGeneticInput


class SetGeneticOCInput(ABCSetInput):
    def __init__(self):
        self.input_parameters = OCGeneticInput()

    def set(self, user_input):
        self.input_parameters.n_chromosomes = int(user_input.genetic.section_dictionary['n_chromosomes'])
        self.input_parameters.n_selected_chr = int(user_input.genetic.section_dictionary['n_selected_chr'])
        self.input_parameters.genetic_algorithm = user_input.genetic.section_dictionary['genetic_algorithm']
        self.input_parameters.amplitude_lim = float(user_input.genetic.section_dictionary['amplitude_lim'])

        self.input_parameters.mate = user_input.mate.section_dictionary['mate']
        self.input_parameters.mate_probability = float(user_input.mate.section_dictionary['mate_probability'])

        self.input_parameters.mutate = user_input.mutate.section_dictionary['mutate']
        self.input_parameters.mutate_probability = float(user_input.mutate.section_dictionary['mutate_probability'])
        self.input_parameters.mutate_mu = float(user_input.mutate.section_dictionary['mutate_mu'])
        self.input_parameters.mutate_starting_sigma = float(user_input.mutate.section_dictionary['mutate_starting_sigma'])
        self.input_parameters.select = user_input.select.section_dictionary['select']

        self.input_parameters.string_file_config = user_input.string_file_config
