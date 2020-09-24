from read_and_set.input.AlphaInput import AlphaInput
from read_and_set.set.ABCSetInput import ABCSetInput


class SetAlphaInput(ABCSetInput):
    def __init__(self):
        super().__init__()
        self.input_parameters = AlphaInput()
        self.read_restart = None

    def set(self, user_input):
        self.input_parameters.dt = float(user_input.sys.section_dictionary['dt'])
        self.input_parameters.nstep = int(user_input.sys.section_dictionary['nstep'])
        self.input_parameters.alpha_type = user_input.oc.section_dictionary['alpha']
        self.input_parameters.read_name = user_input.sys.section_dictionary['folder'] + \
                                                   user_input.oc.section_dictionary['alpha_file']