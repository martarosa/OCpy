import numpy as np

from propagator import PropagatorOCRabitz as rabitzI
from read_and_set.read import auxiliary_functions as af

from field.PropagatorFieldOC import PropagatorFieldOC
from ABCOCIterator import ABCOCIterator, OCIteratorParameters, SimulationParameters



class OCRabitzIterator(ABCOCIterator):
    def __init__(self):
        super().__init__()
        self.par = OCIteratorParameters()
        self.simulation_par = SimulationParameters()

        self.convergence_t = None
        self.J = None
        self.field_psi_matrix = None
        self.psi_coeff_t = None
        self.dict_out = {}

        #Rabitz
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
                for i in range(self.simulation_par.nstep, 0, -1):
                    if i == self.simulation_par.nstep:
                        if self.rabitz_iterator == 'rabitzi':
                            self.prop_chi.propagator_terms.mol.wf.set_wf(af.apply_projection(self.prop_psi.propagator_terms.mol.wf.ci,
                                                                                             self.par.target_state),
                                                                         0)
                        elif self.rabitz_iterator == 'rabitzii':
                            self.prop_chi.propagator_terms.mol.wf.set_wf(self.par.target_state, 0)


                        self.chi_coeff_t[i] = self.prop_chi.propagator_terms.mol.wf.ci

                    self.field_chi_vector_t[i - 1, 1:] = self.prop_field.propagate_field_OC_Rabitz(
                        self.par.psi_coeff_t[i],
                        self.prop_chi.propagator_terms.mol.wf.ci,
                        self.prop_psi.propagator_terms.mol.muT + self.prop_psi.propagator_terms.pcm.par.muLF,
                        self.par.alpha_t[i - 1])
                    self.prop_chi.propagate_one_step(i, self.simulation_par.dt, self.prop_field.field_dt,
                                                     self.par.psi_coeff_t[i])

                    self.chi_coeff_t[i - 1] = self.prop_chi.propagator_terms.mol.wf.ci

                for i in range(self.simulation_par.nstep):
                    if i == 0:
                        self.prop_psi.propagator_terms.mol.wf.set_wf(self.initial_c0, 0)
                        self.par.psi_coeff_t[i] = self.prop_psi.propagator_terms.mol.wf.ci

                    self.field_psi_matrix[i,1:] = self.prop_field.propagate_field_OC_Rabitz(
                        self.prop_psi.propagator_terms.mol.wf.ci,
                        self.chi_coeff_t[i],
                        self.prop_psi.propagator_terms.mol.muT + self.prop_psi.propagator_terms.pcm.par.muLF,
                        self.par.alpha_t[i])
                        # qui invece scorre
                    self.prop_psi.propagate_one_step(i, self.simulation_par.dt, self.prop_field.field_dt)
                    self.par.psi_coeff_t[i + 1] = self.prop_psi.propagator_terms.mol.wf.ci  # coefficients are stored

                self.check_convergence( )

            else:
                self.par.psi_coeff_t = self.prop_psi.propagate_n_step(self.simulation_par.dt,
                                                                      self.simulation_par.nstep,
                                                                      self.field_psi_matrix[:,1:])


    def check_convergence(self):
        J_prev_tmp = np.copy(self.J)
        self.calc_J()
        self.par.convergence_t = self.J - J_prev_tmp


    def calc_J(self):
        if self.rabitz_iterator == "rabitzi":
            self.par.J = np.real(af.projector_mean_value(self.prop_psi.propagator_terms.mol.wf.ci, self.par.target_state) \
                                 - af.alpha_field_J_integral(self.field_psi_matrix, self.par.alpha_t, self.simulation_par.dt))
        elif self.rabitz_iterator == "rabitzii":
            self.par.J = np.real(2 * np.real(np.dot(self.par.target_state, self.prop_psi.propagator_terms.mol.wf.ci) \
                                             - af.alpha_field_J_integral(self.field_psi_matrix, self.par.alpha_t, self.simulation_par.dt)))


    def init_output_dictionary(self):
        self.dict_out['log_file'] = self.get_log_file_out
        self.dict_out['final_pop'] = self.get_final_pop
        self.dict_out['pop_t'] = self.get_pop_t
        self.dict_out['field_t'] = self.get_field_t

    def init(self, oc_input, molecule, starting_field, pcm, alpha_t):
        self.par.target_state = oc_input.target_state
        self.par.alpha_t = alpha_t
        self.J = 99999
        self.convergence_t = 99999

        self.simulation_par.dt = oc_input.dt
        self.simulation_par.nstep = oc_input.nstep

        self.field_psi_matrix = np.copy(starting_field.get_full_field())
        self.psi_coeff_t = np.zeros([self.simulation_par.nstep + 1, molecule.wf.n_ci], dtype=complex)

        self.rabitz_iterator = oc_input.oc_iterator_name

        self.prop_psi.set_propagator(molecule, pcm)
        self.prop_chi.set_propagator(molecule, pcm)
        self.chi_coeff_t = np.zeros([self.simulation_par.nstep + 1, self.prop_psi.propagator_terms.mol.wf.n_ci], dtype=complex)
        self.field_chi_vector_t = np.zeros([self.simulation_par.nstep, 4], dtype=complex)
        self.field_chi_vector_t[:,0] = self.field_psi_matrix[:, 0]

        self.initial_c0 = self.prop_psi.propagator_terms.mol.wf.ci
        self.psi_coeff_t[0] =  self.initial_c0
        self.init_output_dictionary()


# methods inside dictionary return things to be saved
    def get_log_file_out(self, dt):
        integral_field = np.real(af.field_J_integral(self.field_psi_matrix, self.simulation_par.dt))
        norm_proj = np.real(af.projector_mean_value(self.prop_psi.propagator_terms.mol.wf.ci, self.par.target_state)
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





