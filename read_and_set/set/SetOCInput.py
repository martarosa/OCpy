from read_and_set.read import auxiliary_functions as af
from read_and_set.set.ABCSetInput import ABCSetInput
from read_and_set.input.OCInput import OCInput


class SetOCInput(ABCSetInput):
    def __init__(self):
        self.input_parameters = OCInput()

    def set(self, user_input):
        self.input_parameters.oc_iterator_name = user_input.sys.section_dictionary['oc_algorithm']
        self.input_parameters.alpha = user_input.oc.section_dictionary['alpha']
        self.input_parameters.alpha0 = float(user_input.oc.section_dictionary['alpha0'])
        self.input_parameters.n_iterations = int(user_input.oc.section_dictionary['n_iterations'])
        self.input_parameters.convergence_thr = float(user_input.oc.section_dictionary['convergence_thr'])
        
        self.input_parameters.propagation_type = user_input.sys.section_dictionary['propagation_type']
        self.input_parameters.optimization_method = user_input.sys.section_dictionary['scipy_method']
        self.input_parameters.restart = user_input.oc.section_dictionary['restart']

        self.input_parameters.nstep = int(user_input.sys.section_dictionary['nstep'])
        self.input_parameters.dt = float(user_input.sys.section_dictionary['dt'])
        self.input_parameters.target_state = af.normalize_vector([float(i) for i in user_input.oc.section_dictionary["target_state"].split(' ')])
        self.input_parameters.iterator_config_file_name = user_input.oc.section_dictionary['iterator_config_file']





