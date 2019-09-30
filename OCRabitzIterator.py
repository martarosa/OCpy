import numpy as np


import PropagatorOCRabitz as rabitzI
import auxiliary_functions as af


from Field import PropagatorFieldOC
from OCIterator import OCIterator



class OCRabitzIterator(OCIterator):
    def __init__(self):
        super().__init__()
        self.nstep = None
        self.dt = None
        self.target_state = None
        self.alpha_t = None

        self.convergence_t =  99999
        self.J = 99999
        self.initial_c0 = None


        self.rabitz_iterator = None
        self.prop_psi = rabitzI.PropagatorOCfwd()
        self.prop_chi = rabitzI.PropagatorOCbwd()
        self.prop_field = PropagatorFieldOC()

        # storage vector optimal control
        self.psi_coeff_t = None
        self.chi_coeff_t = None
        self.field_psi_matrix = None
        self.field_chi_vector_t = None


        self.dict_out = {}
        self.dict_restart = {}


    def iterate(self, current_iteration):
            if current_iteration != 0:
                for i in range(self.nstep, 0, -1):
                    if i == self.nstep:
                        if self.rabitz_iterator == 'rabitzi':
                            self.prop_chi.propagator_terms.mol.wf.set_wf(af.apply_projection(self.prop_psi.propagator_terms.mol.wf.ci,
                                                                            self.target_state),
                                                        0)
                        elif self.rabitz_iterator == 'rabitzii':
                            self.prop_chi.propagator_terms.mol.wf.set_wf(self.target_state, 0)


                        self.chi_coeff_t[i] = self.prop_chi.propagator_terms.mol.wf.ci
                    self.field_chi_vector_t[i - 1] = self.prop_field.propagate_field_OC_Rabitz(
                        self.psi_coeff_t[i],
                        self.prop_chi.propagator_terms.mol.wf.ci,
                        self.prop_psi.propagator_terms.mol.muT + self.prop_psi.propagator_terms.pcm.muLF,
                        self.alpha_t[i - 1])
                    self.prop_chi.propagate_one_step(self.prop_field.field_dt,
                                                     self.psi_coeff_t[i])

                    self.chi_coeff_t[i - 1] = self.prop_chi.propagator_terms.mol.wf.ci

                for i in range(self.nstep):
                    if i == 0:
                        self.prop_psi.propagator_terms.mol.wf.set_wf(self.initial_c0, 0)
                        self.psi_coeff_t[i] = self.prop_psi.propagator_terms.mol.wf.ci

                    self.field_psi_matrix[i] = self.prop_field.propagate_field_OC_Rabitz(
                        self.prop_psi.propagator_terms.mol.wf.ci,
                        self.chi_coeff_t[i],
                        self.prop_psi.propagator_terms.mol.muT + self.prop_psi.propagator_terms.pcm.muLF,
                        self.alpha_t[i])
                        # qui invece scorre
                    self.prop_psi.propagate_one_step(self.prop_field.field_dt)
                    self.psi_coeff_t[i + 1] = self.prop_psi.propagator_terms.mol.wf.ci  # coefficients are stored

                self.check_convergence( )

            else:
                self.psi_coeff_t = self.prop_psi.propagate_n_step(self.nstep, self.field_psi_matrix)




    def check_convergence(self):
        J_prev_tmp = np.copy(self.J)
        self.calc_J()
        self.convergence_t = self.J - J_prev_tmp

    def calc_J(self):
        if self.rabitz_iterator == "rabitzi":
            self.J = np.real(af.projector_mean_value(self.prop_psi.propagator_terms.mol.wf.ci, self.target_state) \
                             - af.alpha_field_J_integral(self.field_psi_matrix, self.alpha_t, self.dt))
        elif self.rabitz_iterator == "rabitzii":
            self.J = np.real(2 * np.real(np.dot(self.target_state, self.prop_psi.propagator_terms.mol.wf.ci) \
                                         - af.alpha_field_J_integral(self.field_psi_matrix, self.alpha_t, self.dt)))





    def init_output_dictionary(self):
        self.dict_out['log_file'] = self.get_log_file_out
        self.dict_out['final_pop'] = self.get_final_pop
        self.dict_out['pop_t'] = self.get_pop_t
        self.dict_out['field_t'] = self.get_field_t

    def init_oc_iterator(self, molecule, starting_field, pcm, oc_parameters, alpha_t):
        self.nstep = oc_parameters.nstep
        self.dt = oc_parameters.dt
        self.target_state = oc_parameters.target_state
        self.rabitz_iterator = oc_parameters.oc_iterator_name
        self.alpha_t = alpha_t

        self.prop_psi.set_propagator(self.dt, molecule, pcm)
        self.prop_chi.set_propagator(-1 * self.dt, molecule, pcm)

        self.field_psi_matrix = np.copy(starting_field.field)

        self.psi_coeff_t = np.zeros([self.nstep + 1, self.prop_psi.propagator_terms.mol.wf.n_ci], dtype=complex)
        self.chi_coeff_t = np.zeros([self.nstep + 1, self.prop_psi.propagator_terms.mol.wf.n_ci], dtype=complex)
        self.field_chi_vector_t = np.zeros([self.nstep, 3], dtype=complex)

        self.initial_c0 = self.prop_psi.propagator_terms.mol.wf.ci
        self.psi_coeff_t[0] =  self.initial_c0
        self.init_output_dictionary()



# methods inside dictionary return things to be saved
    def get_log_file_out(self):
        integral_field = np.real(af.field_J_integral(self.field_psi_matrix, self.dt))
        norm_proj = np.real(af.projector_mean_value(self.prop_psi.propagator_terms.mol.wf.ci, self.target_state)
                            /(np.dot(self.prop_psi.propagator_terms.mol.wf.ci, np.conj(self.prop_psi.propagator_terms.mol.wf.ci))))
        return np.array([self.convergence_t, self.J, norm_proj, integral_field])

    def get_final_pop(self):
        final_pop = np.real(af.population_from_wf_vector(self.prop_psi.propagator_terms.mol.wf.ci))
        return final_pop

    def get_pop_t(self):
        pop_t = np.real(af.population_from_wf_matrix(self.psi_coeff_t))
        return pop_t

    def get_field_t(self):
        return self.field_psi_matrix

    def get_convergence(self):
        return self.convergence_t

    def get_restart(self):
        return self.get_field_t()





