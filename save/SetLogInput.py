import numpy as np

from SetInput import SetInput
from save.LogHeaderInput import LogHeaderInput


class SetLogInput(SetInput):
    def __init__(self):
        super().__init__()
        self.input_parameters = LogHeaderInput()

    def set(self, user_input):
        self.input_parameters.dt = user_input.sys.section_dictionary['dt']
        self.input_parameters.env = user_input.env.section_dictionary['env']

        self.input_parameters.restart = user_input.oc.section_dictionary['restart']
        self.input_parameters.target_state = user_input.oc.section_dictionary["target_state"]
        self.input_parameters.alpha = user_input.oc.section_dictionary['alpha']
        self.input_parameters.alpha0 = user_input.oc.section_dictionary['alpha0']

        self.input_parameters.field_type = user_input.field.section_dictionary['field_type']
        self.input_parameters.fi = np.array2string(user_input.field.section_dictionary['fi'])
        self.input_parameters.sigma = user_input.field.section_dictionary['sigma']
        self.input_parameters.omega = np.array2string(user_input.field.section_dictionary['omega'])
        self.input_parameters.t0 = user_input.field.section_dictionary['t0']
