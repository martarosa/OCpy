import numpy as np

from propagator import PropagatorOCRabitz as rabitzI
from read import auxiliary_functions as af

from field.Field import PropagatorFieldOC
from OCIterator import OCIterator, OCIteratorAttr



class OCRabitzIterator(OCIterator):
    def __init__(self):
        super().__init__()
        self.class_attributes = OCIteratorAttr()

        self.initial_c0 = None


        self.rabitz_iterator = None
        self.prop_psi = rabitzI.PropagatorOCfwd()
        self.prop_chi = rabitzI.PropagatorOCbwd()
        self.prop_field = PropagatorFieldOC()

        # storage vector optimal control
        self.chi_coeff_t = None
        self.field_chi_vector_t = None




    def iterate(self, current_iteration):
            if current_iteration != 0:
                for i in range(self.class_attributes.nstep, 0, -1):
                    if i == self.class_attributes.nstep:
                        if self.rabitz_iterator == 'rabitzi':
                            self.prop_chi.propagator_terms.mol.wf.set_wf(af.apply_projection(self.prop_psi.propagator_terms.mol.wf.ci,
                                                                            self.class_attributes.target_state),
                                                        0)
                        elif self.rabitz_iterator == 'rabitzii':
                            self.prop_chi.propagator_terms.mol.wf.set_wf(self.class_attributes.target_state, 0)


                        self.chi_coeff_t[i] = self.prop_chi.propagator_terms.mol.wf.ci
                    self.field_chi_vector_t[i - 1] = self.prop_field.propagate_field_OC_Rabitz(
                        self.class_attributes.psi_coeff_t[i],
                        self.prop_chi.propagator_terms.mol.wf.ci,
                        self.prop_psi.propagator_terms.mol.muT + self.prop_psi.propagator_terms.pcm.muLF,
                        self.class_attributes.alpha_t[i - 1])
                    self.prop_chi.propagate_one_step(self.prop_field.field_dt,
                                                     self.class_attributes.psi_coeff_t[i])

                    self.chi_coeff_t[i - 1] = self.prop_chi.propagator_terms.mol.wf.ci

                for i in range(self.class_attributes.nstep):
                    if i == 0:
                        self.prop_psi.propagator_terms.mol.wf.set_wf(self.initial_c0, 0)
                        self.class_attributes.psi_coeff_t[i] = self.prop_psi.propagator_terms.mol.wf.ci

                    self.class_attributes.field_psi_matrix[i] = self.prop_field.propagate_field_OC_Rabitz(
                        self.prop_psi.propagator_terms.mol.wf.ci,
                        self.chi_coeff_t[i],
                        self.prop_psi.propagator_terms.mol.muT + self.prop_psi.propagator_terms.pcm.muLF,
                        self.class_attributes.alpha_t[i])
                        # qui invece scorre
                    self.prop_psi.propagate_one_step(self.prop_field.field_dt)
                    self.class_attributes.psi_coeff_t[i + 1] = self.prop_psi.propagator_terms.mol.wf.ci  # coefficients are stored

                self.check_convergence( )

            else:
                self.class_attributes.psi_coeff_t = self.prop_psi.propagate_n_step(self.class_attributes.nstep, self.class_attributes.field_psi_matrix)






    def check_convergence(self):
        J_prev_tmp = np.copy(self.class_attributes.J)
        self.calc_J()
        self.class_attributes.convergence_t = self.class_attributes.J - J_prev_tmp


    def calc_J(self):
        if self.rabitz_iterator == "rabitzi":
            self.class_attributes.J = np.real(af.projector_mean_value(self.prop_psi.propagator_terms.mol.wf.ci, self.class_attributes.target_state) \
                             - af.alpha_field_J_integral(self.class_attributes.field_psi_matrix, self.class_attributes.alpha_t, self.class_attributes.dt))
        elif self.rabitz_iterator == "rabitzii":
            self.class_attributes.J = np.real(2 * np.real(np.dot(self.class_attributes.target_state, self.prop_psi.propagator_terms.mol.wf.ci) \
                                         - af.alpha_field_J_integral(self.class_attributes.field_psi_matrix, self.class_attributes.alpha_t, self.class_attributes.dt)))





    def init_output_dictionary(self):
        self.class_attributes.dict_out['log_file'] = self.get_log_file_out
        self.class_attributes.dict_out['final_pop'] = self.get_final_pop
        self.class_attributes.dict_out['pop_t'] = self.get_pop_t
        self.class_attributes.dict_out['field_t'] = self.get_field_t

    def init(self, oc_parameters, iterator_parameters, molecule, starting_field, pcm, alpha_t):
        self.class_attributes.nstep = oc_parameters.nstep
        self.class_attributes.dt = oc_parameters.dt
        self.class_attributes.target_state = oc_parameters.target_state
        self.class_attributes.alpha_t = alpha_t
        self.class_attributes.J = 99999
        self.class_attributes.convergence_t = 99999
        self.class_attributes.field_psi_matrix = np.copy(starting_field.field)
        self.class_attributes.psi_coeff_t = np.zeros([self.class_attributes.nstep + 1, molecule.wf.n_ci], dtype=complex)

        self.rabitz_iterator = oc_parameters.oc_iterator_name
        self.prop_psi.set_propagator(self.class_attributes.dt, molecule, pcm)
        self.prop_chi.set_propagator(-1 * self.class_attributes.dt, molecule, pcm)
        self.chi_coeff_t = np.zeros([self.class_attributes.nstep + 1, self.prop_psi.propagator_terms.mol.wf.n_ci], dtype=complex)
        self.field_chi_vector_t = np.zeros([self.class_attributes.nstep, 3], dtype=complex)

        self.initial_c0 = self.prop_psi.propagator_terms.mol.wf.ci
        self.class_attributes.psi_coeff_t[0] =  self.initial_c0
        self.init_output_dictionary()

        # iterator_parameters = rabitz EMPTY




# methods inside dictionary return things to be saved
    def get_log_file_out(self):
        integral_field = np.real(af.field_J_integral(self.class_attributes.field_psi_matrix, self.class_attributes.dt))
        norm_proj = np.real(af.projector_mean_value(self.prop_psi.propagator_terms.mol.wf.ci, self.class_attributes.target_state)
                            /(np.dot(self.prop_psi.propagator_terms.mol.wf.ci, np.conj(self.prop_psi.propagator_terms.mol.wf.ci))))
        return np.array([self.class_attributes.convergence_t, self.class_attributes.J, norm_proj, integral_field])

    def get_final_pop(self):
        final_pop = np.real(af.population_from_wf_vector(self.prop_psi.propagator_terms.mol.wf.ci))
        return final_pop

    def get_pop_t(self):
        pop_t = np.real(af.population_from_wf_matrix(self.class_attributes.psi_coeff_t))
        return pop_t

    def get_field_t(self):
        return self.class_attributes.field_psi_matrix

    def get_convergence(self):
        return self.class_attributes.convergence_t

    def get_restart(self):
        return self.get_field_t()





