import numpy as np

from propagator import PropagatorsEulero as prop
from read_and_set.read import auxiliary_functions as af

from ABCOCIterator import ABCOCIterator, OCIteratorParameters, SimulationParameters




class Eulero1PropagationIterator(ABCOCIterator):
    def __init__(self):
        super().__init__()
        self.par = OCIteratorParameters()
        self.simulation_par = SimulationParameters()

        self.convergence_t = None
        self.J = None
        self.field_psi_matrix = None
        self.psi_coeff_t = None
        self.dict_out = {}

        self.prop_psi = prop.PropagatorEulero1Order()


    def iterate(self, current_iteration):
        self.par.psi_coeff_t = self.prop_psi.propagate_n_step(self.simulation_par.dt,
                                                              self.simulation_par.nstep,
                                                              self.field_psi_matrix)


    def init(self, oc_input, molecule, starting_field, env, alpha_t):
        self.par.target_state = oc_input.target_state
        self.par.alpha_t = alpha_t
        self.J = 99999
        self.convergence_t = 99999

        self.simulation_par.dt = oc_input.dt
        self.simulation_par.nstep = oc_input.nstep

        self.prop_psi.set_propagator(molecule, env)
        self.field_psi_matrix = np.copy(starting_field.get_full_field())

        self.psi_coeff_t = np.zeros([self.simulation_par.nstep + 1, self.prop_psi.propagator_terms.mol.wf.n_ci], dtype=complex)
        self.init_output_dictionary()


    def init_output_dictionary(self):
        self.dict_out['pop_t'] = self.get_pop_t

    def get_pop_t(self):
        pop_t = np.real(af.population_from_wf_matrix(self.psi_coeff_t))
        return pop_t


    def get_restart(self):
        return self.field_psi_matrix


    def check_convergence(self):
        pass

    def calc_J(self):
        pass




class Eulero2PropagationIterator(ABCOCIterator):
    def __init__(self):
        super().__init__()
        self.par = OCIteratorParameters()
        self.simulation_par = SimulationParameters()

        self.convergence_t = None
        self.J = None
        self.field_psi_matrix = None
        self.psi_coeff_t = None
        self.dict_out = {}
        self.prop_psi = prop.PropagatorEulero2Order()


    def iterate(self, current_iteration):
        self.par.psi_coeff_t = self.prop_psi.propagate_n_step(self.simulation_par.dt,
                                                              self.simulation_par.nstep,
                                                              self.field_psi_matrix)

    def init(self, oc_input, molecule, starting_field, env, alpha_t):
        self.par.target_state = oc_input.target_state
        self.par.alpha_t = alpha_t
        self.J = 99999
        self.convergence_t = 99999

        self.simulation_par.dt = oc_input.dt
        self.simulation_par.nstep = oc_input.nstep

        self.prop_psi.set_propagator(molecule, env)
        self.field_psi_matrix = np.copy(starting_field.get_full_field)

        self.psi_coeff_t = np.zeros([self.simulation_par.nstep + 1, self.prop_psi.propagator_terms.mol.wf.n_ci], dtype=complex)
        self.init_output_dictionary()


    def init_output_dictionary(self):
        self.dict_out['pop_t'] = self.get_pop_t

    def get_pop_t(self):
        pop_t = np.real(af.population_from_wf_matrix(self.psi_coeff_t))
        return pop_t

    def get_restart(self):
        return self.field_psi_matrix


    def check_convergence(self):
        pass


    def calc_J(self):
        pass



