import numpy as np

from propagator import PropagatorsEulero as prop
from read import auxiliary_functions as af

from OCIterator import OCIterator




class Eulero1PropagationIterator(OCIterator):
    def __init__(self):
        super().__init__()
        self.nstep = None
        self.dt = None
        self.alpha_t = None
        self.convergence_t =  99999
        self.J = 99999

        self.field_psi_matrix = None
        self.psi_coeff_t = None

        self.prop_psi = prop.PropagatorEulero1Order()
        self.dict_out = {}
        self.dict_restart = {}

    def iterate(self, current_iteration):
        self.psi_coeff_t = self.prop_psi.propagate_n_step(self.nstep, self.field_psi_matrix)


    def init_oc_iterator(self, oc_parameters, iterator_parameters, molecule, starting_field, env, alpha_t):
        self.nstep = oc_parameters.nstep
        self.dt = oc_parameters.dt
        self.alpha_t = alpha_t

        self.prop_psi.set_propagator(self.dt, molecule, env)
        self.field_psi_matrix = np.copy(starting_field.field)

        self.psi_coeff_t = np.zeros([self.nstep + 1, self.prop_psi.propagator_terms.mol.wf.n_ci], dtype=complex)
        self.init_output_dictionary()


    def init_output_dictionary(self):
        self.dict_out['pop_t'] = self.get_pop_t


    def get_pop_t(self):
        pop_t = np.real(af.population_from_wf_matrix(self.psi_coeff_t))
        return pop_t

    def get_restart(self):
        return self.field_psi_matrix


class Eulero2PropagationIterator(OCIterator):
    def __init__(self):
        super().__init__()
        self.nstep = None
        self.dt = None
        self.convergence_t =  99999
        self.J = 99999

        self.field_psi_matrix = None
        self.psi_coeff_t = None

        self.prop_psi = prop.PropagatorEulero2Order()
        self.dict_out = {}
        self.dict_restart = {}

    def iterate(self, current_iteration):
        self.psi_coeff_t = self.prop_psi.propagate_n_step(self.nstep, self.field_psi_matrix)

    def init_oc_iterator(self, oc_parameters, iterator_parameters, molecule, starting_field, env,  alpha_t):
        self.nstep = oc_parameters.nstep
        self.dt = oc_parameters.dt

        self.prop_psi.set_propagator(self.dt, molecule, env)
        self.field_psi_matrix = np.copy(starting_field.field)

        self.psi_coeff_t = np.zeros([self.nstep + 1, self.prop_psi.propagator_terms.mol.wf.n_ci], dtype=complex)
        self.init_output_dictionary()


    def init_output_dictionary(self):
        self.dict_out['pop_t'] = self.get_pop_t

    def get_pop_t(self):
        pop_t = np.real(af.population_from_wf_matrix(self.psi_coeff_t))
        return pop_t

    def get_restart(self):
        return self.field_psi_matrix




