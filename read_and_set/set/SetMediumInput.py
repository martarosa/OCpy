from read_and_set.input.MediumInput import MediumInput
from read_and_set.read.external_output.ReadOutputGaussian import ReadOutputGaussian
from read_and_set.set.ABCSetInput import ABCSetInput


class SetMediumInput(ABCSetInput):
    def __init__(self):
        super().__init__()
        self.input_parameters = MediumInput()

    def set(self, user_input):
        read_output = ReadOutputGaussian()
        self.input_parameters.medium = user_input.medium.section_dictionary['medium']
        if self.input_parameters.medium != 'vac':
            self.input_parameters.cavity = read_output.read_cavity_tesserae(user_input.sys.section_dictionary['folder'] +
                                                                            user_input.medium.section_dictionary['name_file_cavity'])
            if self.input_parameters.medium == 'sol':
                 self.input_parameters.Qnn_reactionfield = read_output.read_Q_matrix(user_input.sys.section_dictionary['folder'] +
                                                                                 user_input.medium.section_dictionary['name_q_reaction_field_dyn'])
                 self.input_parameters.Qnn_localfield = read_output.read_Q_matrix(user_input.sys.section_dictionary['folder'] +
                                                                              user_input.medium.section_dictionary['name_q_local_field_dyn'])