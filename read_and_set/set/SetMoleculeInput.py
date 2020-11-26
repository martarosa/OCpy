from read_and_set.input.MoleculeInput import MoleculeInput
from read_and_set.read.external_output.ReadOutputGaussian import ReadOutputGaussian
from read_and_set.read.external_output.ReadOutputPsi4 import ReadOutputPsi4
from read_and_set.set.ABCSetInput import ABCSetInput
import read_and_set.read.auxiliary_functions as af


class SetMoleculeInput(ABCSetInput):
    def __init__(self):
        super().__init__()
        self.input_parameters = MoleculeInput()


    def set(self, user_input):
        if user_input.oc.section_dictionary['oc_problem'] == 'ground_state':
            read_output = ReadOutputPsi4().read_psi4_output_file(user_input.sys.section_dictionary['folder'] +
                                                               user_input.wf.section_dictionary['name_quantum_chemistry'])
            
#### temporaneo, magari è più carino se faccio scaricare un dizionario e ci sono le keyword #####            
            self.input_parameters.ao_coefficients = read_output[0]
            self.input_parameters.doubly_occupied_mo = read_output[1]
            self.input_parameters.nmo = read_output[2]    
            self.input_parameters.two_electron_int = read_output[5]
            self.input_parameters.ao_kinetic_integral = read_output[3]
            self.input_parameters.ao_potential_integral = read_output[4]
            
#### temporaneo, bisogna che l'init della molecola non si disturbi se non ha questi parametri ###
            read_output = ReadOutputGaussian()
            self.input_parameters.wf_ci = read_output.read_ci0(user_input.sys.section_dictionary['folder'] +
                                                               user_input.wf.section_dictionary['name_ci'])
    
            self.input_parameters.en_ci = read_output.read_en_ci0(user_input.sys.section_dictionary['folder'] +
                                                                  user_input.wf.section_dictionary['name_ei'])
            self.input_parameters.muT = read_output.read_muT(user_input.sys.section_dictionary['folder'] +
                                                             user_input.wf.section_dictionary['name_mut'], self.input_parameters.wf_ci.size)
        
        elif user_input.sys.oc_problem['oc_problem'] == 'optical_excitation':
            read_output = ReadOutputGaussian()
            self.input_parameters.wf_ci = read_output.read_ci0(user_input.sys.section_dictionary['folder'] +
                                                               user_input.wf.section_dictionary['name_ci'])
    
            self.input_parameters.en_ci = read_output.read_en_ci0(user_input.sys.section_dictionary['folder'] +
                                                                  user_input.wf.section_dictionary['name_ei'])
            self.input_parameters.muT = read_output.read_muT(user_input.sys.section_dictionary['folder'] +
                                                             user_input.wf.section_dictionary['name_mut'], self.input_parameters.wf_ci.size)
        else:
            af.exit_error("WARNING: a specific control problem keyword must be included in the input file.")
        if user_input.medium.section_dictionary['medium'] != 'vac':
            self.input_parameters.Vijn = read_output.read_V(user_input.sys.section_dictionary['folder'] +
                                                            user_input.medium.section_dictionary['name_vij'],
                                                            self.input_parameters.wf_ci.size)