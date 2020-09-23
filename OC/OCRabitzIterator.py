import numpy as np
import pandas as pd
from copy import deepcopy

import dictionaries.PropagatorDictionaries
from propagator import PropagatorOCRabitz as prop
from read_and_set.read import auxiliary_functions as af
from dictionaries import SaveDictionaries as dict
from field.PropagatorFieldOC import PropagatorFieldOC
from OC.ABCOCIterator import ABCOCIterator
from parameters.OCIteratorParameters import OCIteratorParameters
from parameters.RabitzParameters import RabitzParameters
from SystemObj import DiscreteTimePar

from field.Field import Func_tMatrix


class OCRabitzIterator(ABCOCIterator):
    def __init__(self):
        super().__init__()
        self.par = OCIteratorParameters()
        self.rabitz_par = RabitzParameters()

        self.discrete_t_par = DiscreteTimePar()


        self.field_psi_matrix = Func_tMatrix()
        self.psi_coeff_t_matrix = Func_tMatrix()

        self.dict_out = {}



        #Rabitz
        self.prop_psi = None

        self.initial_c0 = None
        self.rabitz_iterator = None
        self.prop_chi = prop.PropagatorOCbwd()
        self.prop_field = PropagatorFieldOC()
        # storage vector optimal control
        self.chi_coeff_t_matrix = Func_tMatrix()
        self.field_chi_matrix = Func_tMatrix()



    def iterate(self, current_iteration):
        mut_rabitz_field_prop = deepcopy(self.prop_psi.mol.par.muT)
        if self.prop_psi.medium.par.medium == "sol":
            mut_rabitz_field_prop += self.prop_psi.medium.muLF
        if current_iteration != 0:
            for i in range(self.discrete_t_par.nstep, 0, -1):
                if i == self.discrete_t_par.nstep:
                    if self.rabitz_iterator == 'rabitzi':
                        self.prop_chi.mol.wf.set_wf(af.apply_projection(self.prop_psi.mol.wf.ci,
                                                                        self.par.target_state),
                                                                     0)
                    elif self.rabitz_iterator == 'rabitzii':
                        self.prop_chi.mol.wf.set_wf(self.par.target_state, 0)


                    self.chi_coeff_t_matrix.f_xyz[i] = self.prop_chi.mol.wf.ci


                self.field_chi_matrix.f_xyz[i - 1] = self.prop_field.propagate_field_OC_Rabitz(
                    self.psi_coeff_t_matrix.f_xyz[i],
                    self.prop_chi.mol.wf.ci,
                    mut_rabitz_field_prop,
                    self.par.alpha_t[i - 1])

                self.prop_chi.propagate_one_step(self.discrete_t_par.dt, self.prop_field.field_dt_vector,
                                                 self.psi_coeff_t_matrix.f_xyz[i])

                self.chi_coeff_t_matrix.f_xyz[i - 1] = self.prop_chi.mol.wf.ci

            for i in range(self.discrete_t_par.nstep):
                if i == 0:
                    self.prop_psi.mol.wf.set_wf(self.initial_c0, 0)
                    self.psi_coeff_t_matrix.f_xyz[i] = self.prop_psi.mol.wf.ci

                self.field_psi_matrix.f_xyz[i] = self.prop_field.propagate_field_OC_Rabitz(
                    self.prop_psi.mol.wf.ci,
                    self.chi_coeff_t_matrix.f_xyz[i],
                    mut_rabitz_field_prop,
                    self.par.alpha_t[i])
                self.prop_psi.propagate_one_step(self.discrete_t_par.dt, self.prop_field.field_dt_vector)
                self.psi_coeff_t_matrix.f_xyz[i + 1] = self.prop_psi.mol.wf.ci  # coefficients are stored
            self.check_convergence( )

        else:
            self.psi_coeff_t_matrix = self.prop_psi.propagate_n_step(self.discrete_t_par,
                                                                     self.field_psi_matrix)

            self.chi_coeff_t_matrix = deepcopy(self.psi_coeff_t_matrix)




    def check_convergence(self):
        J_prev_tmp = np.copy(self.par.J)
        self.calc_J()
        self.par.convergence_t = self.par.J - J_prev_tmp

    def calc_J(self):
        if self.rabitz_iterator == "rabitzi":
            self.par.J = np.real(af.projector_mean_value(self.prop_psi.mol.wf.ci, self.par.target_state) \
                             - self.alpha_field_J_integral())
        elif self.rabitz_iterator == "rabitzii":
            self.par.J = np.real(2 * np.real(np.dot(self.par.target_state, self.prop_psi.mol.wf.ci) \
                                         - self.alpha_field_J_integral()))

    def init(self, molecule, starting_field, medium, alpha_t, oc_input, oc_conf):
        self.par.propagator = oc_input.propagator
        self.prop_psi = dictionaries.PropagatorDictionaries.PropagatorDict[oc_input.propagator]()
        if not isinstance(self.prop_psi, prop.PropagatorOCfwd):
            af.exit_error("Error in the initialization of Eulero 1 order propagator")

        self.discrete_t_par.dt = oc_input.dt
        self.discrete_t_par.nstep = oc_input.nstep

        self.par.target_state = oc_input.target_state
        self.rabitz_par.alpha0 = oc_conf.alpha0
        self.par.alpha_t = self.rabitz_par.alpha0 * np.array(alpha_t)
        self.par.J = 99999
        self.par.convergence_t = 99999

        self.field_psi_matrix = deepcopy(starting_field.field)

        self.init_rabitz(oc_input, molecule, medium)

        self.init_output_dictionary()

    # init specific iterator
    def init_rabitz(self, oc_input, molecule, medium):
        self.rabitz_iterator = oc_input.oc_iterator_name
        self.prop_psi.set_propagator(molecule, medium)
        self.prop_chi.set_propagator(molecule, medium)
        self.field_chi_matrix = deepcopy(self.field_psi_matrix)
        self.initial_c0 = self.prop_psi.mol.wf.ci


    def init_output_dictionary(self):
        self.dict_out['log_file'] = self.get_log_file_out
        self.dict_out['final_pop'] = self.get_final_pop
        self.dict_out['pop_t'] = self.get_pop_t
        self.dict_out['field_t'] = self.get_field_t



# methods inside dictionary return things to be saved
    def get_log_file_out(self):
        integral_field = np.real(self.field_J_integral())
        norm_proj = np.real(af.projector_mean_value(self.prop_psi.mol.wf.ci, self.par.target_state)
                            /(np.dot(self.prop_psi.mol.wf.ci, np.conj(self.prop_psi.mol.wf.ci))))
        return np.array([self.par.convergence_t, self.par.J, norm_proj, integral_field])

    def get_final_pop(self):
        final_pop = np.real(af.population_from_wf_vector(self.prop_psi.mol.wf.ci))
        return final_pop

    def get_pop_t(self):
        psi_coeff_t_matrix = np.insert(self.psi_coeff_t_matrix.f_xyz, 0, self.psi_coeff_t_matrix.time_axis, axis = 1)
        pop_t_matrix = np.real(af.population_from_wf_matrix(psi_coeff_t_matrix))
        return pop_t_matrix


    def get_field_t(self):
        field_t_matrix = np.insert(self.field_psi_matrix.f_xyz, 0, self.field_psi_matrix.time_axis, axis = 1)
        return field_t_matrix

    def get_restart(self):
        header = np.array([["#time","field x","field y", "field z"],["###","###","###","###"]])
        out = pd.DataFrame(header)
        field_t_matrix = pd.DataFrame(self.get_field_t())
        out = out.append(field_t_matrix)
        return out



