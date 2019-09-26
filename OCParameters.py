import numpy as np

import PropagatorsEulero as prop
import auxiliary_functions as af


from OCIterator import OCIterator
from Field import Field
from FieldParameters import FieldParameters

class OCParameters:
    def __init__(self):
        self.oc_iterator_name = None
        self.alpha = None
        self.alpha0 = None
        self.n_iterations = None
        self.convergence_thr = None

        self.restart = None

        self.nstep = None
        self.dt = None
        self.target_state = None

        self.iterator_parameters = None





class GeneticParameters:
    def __init__(self):
        self.delta_omega = None
        self.n_chromosomes = None
        self.genes_amplitude_limits = []




class InitGeneticPar():
    def __init__(self):
        self.genetic_parameters = GeneticParameters()

    def init(self, user_input):
        self.genetic_parameters.genes_amplitude_limits.append(user_input.chr.par['amplitude_low_limit'])
        self.genetic_parameters.genes_amplitude_limits.append(user_input.chr.par['amplitude_hight_limit'])
        self.genetic_parameters.delta_omega = user_input.chr.par['delta_omega']
        self.genetic_parameters.n_chromosomes = user_input.chr.par['n_chromosomes']



class Chromosome():
    def __init__(self):
        self.n_amplitude_genes = None  # n_ampitudes
        self.amplitude_genes_vector_cos = None  # amplitudes
        self.amplitude_genes_vector_sin = None  # amplitudes

        self.field = Field()
        self.prop_psi = prop.PropagatorEulero2Order()


    def init_chromosome(self, oc_parameters, molecule, starting_field, pcm):
        self.field = starting_field
        self.n_amplitude_genes = self.field.parameters['fi'].size()
        self.amplitude_genes_vector_sin = np.random.uniform(low=oc_parameters.iterator_parameters.genes_amplitude_limits[0],
                                                        high=oc_parameters.iterator_parameters.genes_amplitude_limits[1],
                                                        size=(self.n_amplitude_genes))
        self.amplitude_genes_vector_cos = np.random.uniform(low=oc_parameters.iterator_parameters.genes_amplitude_limits[0],
                                                            high=oc_parameters.iterator_parameters.genes_amplitude_limits[1],
                                                            size=(self.n_amplitude_genes))

        self.field.parameters['fi'] = self.amplitude_genes_vector_sin
        self.field.parameters['fi_cos'] = self.amplitude_genes_vector_cos
        self.field.chose_field(self.field.field_type)
        self.prop_psi.set_propagator(oc_parameters.dt, molecule, pcm)




class OCGeneticIterator(OCIterator):
    def __init__(self):
        super().__init__()
        self.nstep = None
        self.dt = None
        self.target_state = None
        self.alpha_t = None
        self.convergence_t = 99999
        self.J = 99999

        self.field_psi_matrix = None
        self.psi_coeff_t = None

        self.dict_out = {}
        self.dict_restart = {}

        # genetic

        self.n_chromosomes = None
        self.n_amplitudes_omega_genes = None
        self.genes_amplitude_limits = []
        self.omegas_vector = None
        self.delta = None

        self.chromosomes = []
        self.Js = []


    def iterate(self, current_iteration):
        for i in range(self.n_chromosomes):
            self.chromosomes[i].prop_psi.propagate_n_step(self.nstep, self.chromosomes[i].field.field)
        self.calc_J()
        self.check_convergence()



    def init_oc_iterator(self, molecule, starting_field, pcm, oc_parameters, alpha_t):
        self.nstep = oc_parameters.nstep
        self.dt = oc_parameters.dt
        self.alpha_t = alpha_t
        self.target_state = oc_parameters.target_state

        self.delta = oc_parameters.iterator_parameters.delta
        self.n_chromosomes = oc_parameters.iterator_parameters.n_chromosomes
        self.genes_amplitude_limits = oc_parameters.iterator_parameters.genes_amplitude_limits


     #   self.set_omega_vector(molecule)

        for i in range(self.n_chromosomes):
            tmp = Chromosome()
            tmp.init_chromosome(oc_parameters, molecule, starting_field, pcm)
            self.chromosomes.append(tmp)

#    def init_output_dictionary(self):
#        self.dict_out['log_file'] = self.get_log_file_out
#        self.dict_out['final_pop'] = self.get_final_pop
#        self.dict_out['pop_t'] = self.get_pop_t
#        self.dict_out['field_t'] = self.get_field_t

#    def get_log_file_out(self):
#        integral_field = np.real(af.field_J_integral(self.field_psi_matrix, self.dt))
#        norm_proj = np.real(af.projector_mean_value(self.prop_psi.mol.wf.ci, self.target_state)/(np.dot(self.prop_psi.mol.wf.ci, np.conj(self.prop_psi.mol.wf.ci))))
#        return np.array([self.convergence_t, self.J, norm_proj, integral_field])

#    def get_final_pop(self):
#        final_pop = np.real(af.population_from_wf_vector(self.prop_psi.mol.wf.ci))
#        return final_pop

#    def get_pop_t(self):
#        pop_t = np.real(af.population_from_wf_matrix(self.psi_coeff_t))
#        return pop_t

#    def get_field_t(self):
#        return self.field_psi_matrix

#    def get_restart(self):
#        return self.get_field_t()






    #def set_omega_vector(self, molecule):
    #    self.n_amplitudes_omega_genes = int(molecule.en_ci[-1]) / self.delta + 1
    #    self.omegas_vector = np.zeros(self.n_amplitudes_omega_genes)
    #    self.omegas_vector[0] = 0
    #    for i in range(self.n_amplitudes_omega_genes):
    #        self.omegas_vector[i + 1] = self.omegas_vector[i] + self.delta


    def calc_J(self):
        self.Js = []
        for i in range(self.n_chromosomes):
            self.Js.append(np.real(af.projector_mean_value(self.chromosomes[i].prop_psi.mol.wf.ci, self.target_state) \
                         - af.alpha_field_J_integral(self.chromosomes[i].field.field, self.alpha_t, self.dt)))
        self.Js, self.chromosomes = (list(t) for t in zip(*sorted(zip(self.Js, self.chromosomes), reverse = True)))

        self.J = self.J[0]


    def check_convergence(self):
        J_prev_tmp = np.copy(self.J)
        self.calc_J()
        self.convergence_t = self.J - J_prev_tmp
















    def init_output_dictionary(self):
        self.dict_out['pop_t'] = self.get_pop_t

    def get_pop_t(self):
        pop_t = np.real(af.population_from_wf_matrix(self.psi_coeff_t))
        return pop_t

    def get_restart(self):
        return self.field_psi_matrix
