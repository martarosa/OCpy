from read_and_set.set.ABCSetInput import ABCSetInput
from read_and_set.input.PCMInput import PCMInput
from read_and_set.read.ReadOutputGaussian import ReadOutputGaussian


class SetPCMInput(ABCSetInput):
    def __init__(self):
        super().__init__()
        self.input_parameters = PCMInput()

    def set(self, user_input):
        read_output = ReadOutputGaussian()
        self.input_parameters.env = user_input.env.section_dictionary['env']
        if self.input_parameters.env != 'vac':
            self.input_parameters.cavity = read_output.read_cavity_tesserae(user_input.sys.section_dictionary['folder'] +
                                                                            user_input.env.section_dictionary['name_file_cavity'])
            if self.input_parameters.env == 'sol':
                 self.input_parameters.Qnn_reactionfield = read_output.read_Q_matrix(user_input.sys.section_dictionary['folder'] +
                                                                                 user_input.env.section_dictionary['name_q_tdplas'])
                 self.input_parameters.Qnn_localfield = read_output.read_Q_matrix(user_input.sys.section_dictionary['folder'] +
                                                                              user_input.env.section_dictionary['name_q_local_field'])