import numpy as np

from propagator import PropagatorsEulero as prop
from read import auxiliary_functions as af

from OCIterator import OCIterator, OCIteratorParameters




class Eulero1PropagationIterator(OCIterator):
    def __init__(self):
        super().__init__()
        self.oc_iterator_parameters = OCIteratorParameters()
        self.prop_psi = prop.PropagatorEulero1Order()
#        self.dict_restart = {}

    def iterate(self, current_iteration):
        self.oc_iterator_parameters.psi_coeff_t = self.prop_psi.propagate_n_step(self.oc_iterator_parameters.nstep, self.oc_iterator_parameters.field_psi_matrix)


    def init(self, oc_iterator_parameters, molecule, starting_field, env, alpha_t):
        self.oc_iterator_parameters.nstep = oc_iterator_parameters.nstep
        self.oc_iterator_parameters.dt = oc_iterator_parameters.dt
        self.oc_iterator_parameters.alpha_t = alpha_t
        self.oc_iterator_parameters.convergence_t = 99999
        self.oc_iterator_parameters.J = 99999

        self.prop_psi.set_propagator(self.oc_iterator_parameters.dt, molecule, env)
        self.field_psi_matrix = np.copy(starting_field.field)

        self.psi_coeff_t = np.zeros([self.oc_iterator_parameters.nstep + 1, self.prop_psi.propagator_terms.mol.wf.n_ci], dtype=complex)
        self.init_output_dictionary()


    def init_output_dictionary(self):
        self.oc_iterator_parameters.dict_out['pop_t'] = self.get_pop_t


    def get_pop_t(self):
        pop_t = np.real(af.population_from_wf_matrix(self.psi_coeff_t))
        return pop_t

    def get_restart(self):
        return self.oc_iterator_parameters.field_psi_matrix


class Eulero2PropagationIterator(OCIterator):
    def __init__(self):
        super().__init__()
        self.oc_iterator_parameters = OCIteratorParameters()
        self.prop_psi = prop.PropagatorEulero2Order()


    def iterate(self, current_iteration):
        self.oc_iterator_parameters.psi_coeff_t = self.prop_psi.propagate_n_step(self.oc_iterator_parameters.nstep, self.oc_iterator_parameters.field_psi_matrix)

    def init(self, oc_iterator_parameters, molecule, starting_field, env, alpha_t):
        self.oc_iterator_parameters.nstep = oc_iterator_parameters.nstep
        self.oc_iterator_parameters.dt = oc_iterator_parameters.dt
        self.oc_iterator_parameters.J = 99999
        self.oc_iterator_parameters.convergence_t = 99999

        self.prop_psi.set_propagator(self.oc_iterator_parameters.dt, molecule, env)
        self.oc_iterator_parameters.field_psi_matrix = np.copy(starting_field.field)

        self.oc_iterator_parameters.psi_coeff_t = np.zeros([self.oc_iterator_parameters.nstep + 1, self.prop_psi.propagator_terms.mol.wf.n_ci], dtype=complex)
        self.init_output_dictionary()


    def init_output_dictionary(self):
        self.oc_iterator_parameters.dict_out['pop_t'] = self.get_pop_t

    def get_pop_t(self):
        pop_t = np.real(af.population_from_wf_matrix(self.oc_iterator_parameters.psi_coeff_t))
        return pop_t

    def get_restart(self):
        return self.oc_iterator_parameters.field_psi_matrix




