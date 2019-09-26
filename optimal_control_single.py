import numpy as np
import sys


import PropagatorsEulero as prop
import auxiliary_functions as af
from optimal_control import OC


class SinglePropagation(OC):
    def __init__(self):
        super().__init__()



    def initialize_oc(self, input_par_file, mol, field, env):
        self.initialize_common_oc_parameters(input_par_file)
        self.initialize_propagation_attributes(mol, field, env)
        input_par_file.get_oc_par()["restart"] = "false"
        self.save.creation_save_files(
                                 ["_pop_t.dat"],
                                 ["#Iteration populations \n #fields: n_iteration, time, states population \n"],
                                 input_par_file)


    def initialize_propagation_attributes(self, molecule, field, env):
        self.prop_psi.set_propagator(self.dt, molecule, env)
        self.set_starting_field(field)
        self.psi_coeff_t = np.zeros([self.nstep + 1, self.prop_psi.mol.wf.n_ci], dtype=complex)
        self.field_chi_vector_t = np.zeros([self.nstep, 3], dtype=complex)
        self.set_initial_c0()

    def iterate_oc(self):
        self.psi_coeff_t = self.prop_psi.propagate_n_step(self.nstep,
                                                          self.field_vector_t)
        self.save.save_every_n_iterations(0,
                                          ["_pop_t.dat"],
                                          np.real(self.calc_every_iteration_output()),
                                          "_field_bkp.dat",
                                          self.field_vector_t)

    def calc_J(self):
        pass

    def check_convergence(self):
        pass

    def calc_every_iteration_output(self):
        pop = af.population_from_wf_vector(self.psi_coeff_t[-1])
        return [np.real(pop)]

    def set_initial_c0(self):
        self.initial_c0 = self.prop_psi.mol.wf.Ci()
        if self.state_to_project_on.size != self.initial_c0.size:
            sys.exit(
                "ERROR: Wrong dimension of target state: number of coefficients is different from wavefunction. ")


    def set_starting_field(self, starting_field):
        self.field_vector_t = np.copy(starting_field.get_field())






class SinglePropagation1Eulero(SinglePropagation):
    def __init__(self):
        super().__init__()
        self.prop_psi = prop.PropagatorEulero1Order()






class SinglePropagation2Eulero(SinglePropagation):
    def __init__(self):
        super().__init__()
        self.prop_psi = prop.PropagatorEulero2Order()