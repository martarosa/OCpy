import numpy as np

from propagator import PropagatorsEulero as prop
from read import auxiliary_functions as af

from OCIterator import OCIterator, OCIteratorAttr




class Eulero1PropagationIterator(OCIterator):
    def __init__(self):
        super().__init__()
        self.class_attributes = OCIteratorAttr()
        self.prop_psi = prop.PropagatorEulero1Order()
#        self.dict_restart = {}

    def iterate(self, current_iteration):
        self.class_attributes.psi_coeff_t = self.prop_psi.propagate_n_step(self.class_attributes.nstep, self.class_attributes.field_psi_matrix)


    def init(self, oc_parameters, iterator_parameters, molecule, starting_field, env, alpha_t):
        self.class_attributes.nstep = oc_parameters.nstep
        self.class_attributes.dt = oc_parameters.dt
        self.class_attributes.alpha_t = alpha_t
        self.class_attributes.convergence_t = 99999
        self.class_attributes.J = 99999

        self.prop_psi.set_propagator(self.class_attributes.dt, molecule, env)
        self.field_psi_matrix = np.copy(starting_field.field)

        self.psi_coeff_t = np.zeros([self.class_attributes.nstep + 1, self.prop_psi.propagator_terms.mol.wf.n_ci], dtype=complex)
        self.init_output_dictionary()


    def init_output_dictionary(self):
        self.class_attributes.dict_out['pop_t'] = self.get_pop_t


    def get_pop_t(self):
        pop_t = np.real(af.population_from_wf_matrix(self.psi_coeff_t))
        return pop_t

    def get_restart(self):
        return self.class_attributes.field_psi_matrix


class Eulero2PropagationIterator(OCIterator):
    def __init__(self):
        super().__init__()
        self.class_attributes = OCIteratorAttr()
        self.prop_psi = prop.PropagatorEulero2Order()


    def iterate(self, current_iteration):
        self.class_attributes.psi_coeff_t = self.prop_psi.propagate_n_step(self.class_attributes.nstep, self.class_attributes.field_psi_matrix)

    def init(self, oc_parameters, iterator_parameters, molecule, starting_field, env, alpha_t):
        self.class_attributes.nstep = oc_parameters.nstep
        self.class_attributes.dt = oc_parameters.dt
        self.class_attributes.J = 99999
        self.class_attributes.convergence_t = 99999

        self.prop_psi.set_propagator(self.class_attributes.dt, molecule, env)
        self.class_attributes.field_psi_matrix = np.copy(starting_field.field)

        self.class_attributes.psi_coeff_t = np.zeros([self.class_attributes.nstep + 1, self.prop_psi.propagator_terms.mol.wf.n_ci], dtype=complex)
        self.init_output_dictionary()


    def init_output_dictionary(self):
        self.class_attributes.dict_out['pop_t'] = self.get_pop_t

    def get_pop_t(self):
        pop_t = np.real(af.population_from_wf_matrix(self.class_attributes.psi_coeff_t))
        return pop_t

    def get_restart(self):
        return self.class_attributes.field_psi_matrix




