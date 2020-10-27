import numpy as np
import pandas as pd
from copy import deepcopy

import dictionaries.PropagatorDictionaries as pdict
from read_and_set.read import auxiliary_functions as af


from OC.ABCOCIterator import ABCOCIterator
from parameters.OCIteratorParameters import OCIteratorParameters
from SystemObj import DiscreteTimePar, Func_tMatrix



class SimplePropagationIterator(ABCOCIterator):
    def __init__(self):
        super().__init__()
        self.par = OCIteratorParameters()
        self.discrete_t_par = DiscreteTimePar()

        self.field_psi_matrix = Func_tMatrix()
        self.psi_coeff_t_matrix = Func_tMatrix()

        self.dict_out = {}

        self.prop_psi = None

    def iterate(self, current_iteration):
        self.prop_psi.propagate_n_step(self.discrete_t_par,
                                       self.field_psi_matrix)
        self.psi_coeff_t_matrix = self.prop_psi.wf_matrix_out

    def check_convergence(self):
        pass

    def calc_J(self):
        pass

    def init(self, molecule, starting_field, medium, alpha_t, oc_input, oc_conf = None, prop_conf = None):
        self.par.propagator = oc_input.propagator
        self.prop_psi = pdict.PropagatorDict[oc_input.propagator]()
        self.discrete_t_par.dt = oc_input.dt
        self.discrete_t_par.nstep = oc_input.nstep
        self.par.target_state = oc_input.target_state
        self.par.alpha_t = np.array(alpha_t)
        self.par.J = 99999
        self.par.convergence_t = 99999

        self.prop_psi.set_propagator(molecule, medium, prop_conf)
        self.field_psi_matrix = deepcopy(starting_field.field)

        self.init_output_dictionary()


    def init_output_dictionary(self):
        self.dict_out['pop_t'] = self.get_pop_t

    def get_pop_t(self):
        psi_coeff_t_matrix = np.insert(self.prop_psi.wf_matrix_out.f_xyz, 0, self.prop_psi.wf_matrix_out.time_axis, axis = 1)
        pop_t = np.real(af.population_from_wf_matrix(psi_coeff_t_matrix))
        return pop_t

    def get_restart(self):
        field_t_matrix = np.insert(self.field_psi_matrix.f_xyz, 0, self.field_psi_matrix.time_axis, axis = 1)
        return pd.DataFrame(field_t_matrix)











