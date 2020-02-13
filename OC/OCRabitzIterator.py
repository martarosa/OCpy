import numpy as np
from copy import deepcopy

from propagator import PropagatorOCRabitz as rabitzI
from read_and_set.read import auxiliary_functions as af

from field.PropagatorFieldOC import PropagatorFieldOC
from OC.ABCOCIterator import ABCOCIterator, OCIteratorParameters
from SystemObj import DiscreteTimePar

from field.Field import Func_tMatrix


class OCRabitzIterator(ABCOCIterator):
    def __init__(self):
        super().__init__()
        self.par = OCIteratorParameters()
        self.discrete_t_par = DiscreteTimePar()


        self.field_psi_matrix = Func_tMatrix()
        self.psi_coeff_t_matrix = Func_tMatrix()

        self.dict_out = {}

        #Rabitz
        self.prop_psi = rabitzI.PropagatorOCfwd()

        self.initial_c0 = None
        self.rabitz_iterator = None
        self.prop_chi = rabitzI.PropagatorOCbwd()
        self.prop_field = PropagatorFieldOC()
        # storage vector optimal control
        self.chi_coeff_t_matrix = Func_tMatrix()
        self.field_chi_matrix = Func_tMatrix()


# non va iterate

    def iterate(self, current_iteration):
            if current_iteration != 0:
                for i in range(self.discrete_t_par.nstep, 0, -1):
                    if i == self.discrete_t_par.nstep:
                        if self.rabitz_iterator == 'rabitzi':
                            self.prop_chi.propagator_terms.mol.wf.set_wf(af.apply_projection(self.prop_psi.propagator_terms.mol.wf.ci,
                                                                                             self.par.target_state),
                                                                         0)
                        elif self.rabitz_iterator == 'rabitzii':
                            self.prop_chi.propagator_terms.mol.wf.set_wf(self.par.target_state, 0)


                        self.chi_coeff_t_matrix.f_xyz[i] = self.prop_chi.propagator_terms.mol.wf.ci


                    self.field_chi_matrix.f_xyz[i - 1] = self.prop_field.propagate_field_OC_Rabitz(
                        self.psi_coeff_t_matrix.f_xyz[i],
                        self.prop_chi.propagator_terms.mol.wf.ci,
                        self.prop_psi.propagator_terms.mol.par.muT + self.prop_psi.propagator_terms.pcm.par.muLF,
                        self.par.alpha_t[i - 1])

                    self.prop_chi.propagate_one_step(i, self.discrete_t_par.dt, self.prop_field.field_dt_vector,
                                                     self.psi_coeff_t_matrix.f_xyz[i])

                    self.chi_coeff_t_matrix.f_xyz[i - 1] = self.prop_chi.propagator_terms.mol.wf.ci

                for i in range(self.discrete_t_par.nstep):
                    if i == 0:
                        self.prop_psi.propagator_terms.mol.wf.set_wf(self.initial_c0, 0)
                        self.psi_coeff_t_matrix.f_xyz[i] = self.prop_psi.propagator_terms.mol.wf.ci

                    self.field_psi_matrix.f_xyz[i] = self.prop_field.propagate_field_OC_Rabitz(
                        self.prop_psi.propagator_terms.mol.wf.ci,
                        self.chi_coeff_t_matrix.f_xyz[i],
                        self.prop_psi.propagator_terms.mol.par.muT + self.prop_psi.propagator_terms.pcm.par.muLF,
                        self.par.alpha_t[i])

                    self.prop_psi.propagate_one_step(i, self.discrete_t_par.dt, self.prop_field.field_dt_vector)
                    print("one step")
                    print(self.prop_psi.propagator_terms.mol.wf.ci)
                    self.psi_coeff_t_matrix.f_xyz[i + 1] = self.prop_psi.propagator_terms.mol.wf.ci  # coefficients are stored
                print(self.psi_coeff_t_matrix.f_xyz)
                self.check_convergence( )

            else:
                self.psi_coeff_t_matrix = self.prop_psi.propagate_n_step(self.discrete_t_par,
                                                                         self.field_psi_matrix)

                self.chi_coeff_t_matrix = self.psi_coeff_t_matrix



    def check_convergence(self):
        J_prev_tmp = np.copy(self.par.J)
        self.calc_J()
        self.par.convergence_t = self.par.J - J_prev_tmp

    def calc_J(self):
        if self.rabitz_iterator == "rabitzi":
            self.par.J = np.real(af.projector_mean_value(self.prop_psi.propagator_terms.mol.wf.ci, self.par.target_state) \
                             - self.alpha_field_J_integral())
        elif self.rabitz_iterator == "rabitzii":
            self.par.J = np.real(2 * np.real(np.dot(self.par.target_state, self.prop_psi.propagator_terms.mol.wf.ci) \
                                         - self.alpha_field_J_integral()))

    def init(self, molecule, starting_field, pcm, alpha_t, oc_input, iterator_config_input = None):
        self.discrete_t_par.dt = oc_input.dt
        self.discrete_t_par.nstep = oc_input.nstep

        self.par.target_state = oc_input.target_state
        self.par.alpha_t = alpha_t
        self.par.J = 99999
        self.par.convergence_t = 99999

        self.field_psi_matrix = deepcopy(starting_field.field)

        self.init_rabitz(oc_input, molecule, pcm)

        #self.psi_coeff_t_matrix = np.zeros([self.discrete_t_par.nstep + 1, self.prop_psi.propagator_terms.mol.wf.n_ci + 1], dtype=complex)
        #self.psi_coeff_t_matrix[:, 0] = np.linspace(0, self.discrete_t_par.dt * self.discrete_t_par.nstep, self.discrete_t_par.nstep + 1)
        #self.psi_coeff_t_matrix[0, 1:] =  self.initial_c0

        self.init_output_dictionary()


    def init_output_dictionary(self):
        self.dict_out['log_file'] = self.get_log_file_out
        self.dict_out['final_pop'] = self.get_final_pop
        self.dict_out['pop_t'] = self.get_pop_t
        self.dict_out['field_t'] = self.get_field_t


#init specific iterator

    def init_rabitz(self, oc_input, molecule, pcm):
        self.rabitz_iterator = oc_input.oc_iterator_name
        self.prop_psi.set_propagator(molecule, pcm)
        self.prop_chi.set_propagator(molecule, pcm)
        #self.chi_coeff_t_matrix = np.zeros([self.discrete_t_par.nstep + 1, self.prop_psi.propagator_terms.mol.wf.n_ci + 1], dtype=complex)
        #self.chi_coeff_t_matrix[:, 0] = np.linspace(0, self.discrete_t_par.dt * self.discrete_t_par.nstep, self.discrete_t_par.nstep + 1)
        #self.field_chi_matrix.f_xyz = np.zeros([self.discrete_t_par.nstep, 3], dtype=complex)
        self.field_chi_matrix = self.field_psi_matrix

        self.initial_c0 = self.prop_psi.propagator_terms.mol.wf.ci


# methods inside dictionary return things to be saved
    def get_log_file_out(self):
        integral_field = np.real(self.field_J_integral())
        norm_proj = np.real(af.projector_mean_value(self.prop_psi.propagator_terms.mol.wf.ci, self.par.target_state)
                            /(np.dot(self.prop_psi.propagator_terms.mol.wf.ci, np.conj(self.prop_psi.propagator_terms.mol.wf.ci))))
        return np.array([self.par.convergence_t, self.par.J, norm_proj, integral_field])

    def get_final_pop(self):
        final_pop = np.real(af.population_from_wf_vector(self.prop_psi.propagator_terms.mol.wf.ci))
        return final_pop

    def get_pop_t(self):
        psi_coeff_t_matrix = np.insert(self.psi_coeff_t_matrix.f_xyz, 0, self.psi_coeff_t_matrix.time_axis, axis = 1)
        pop_t_matrix = np.real(af.population_from_wf_matrix(psi_coeff_t_matrix))
        return pop_t_matrix


    def get_field_t(self):
        field_t_matrix = np.insert(self.field_psi_matrix.f_xyz, 0, self.field_psi_matrix.time_axis, axis = 1)
        return field_t_matrix

    def get_restart(self):
        return self.get_field_t()



