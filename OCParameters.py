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
        self.n_chromosomes = None
        self.genes_amplitude_limits = []



class InitGeneticPar():
    def __init__(self):
        self.genetic_parameters = GeneticParameters()

    def init(self, user_input):
        self.genetic_parameters.genes_amplitude_limits.append(user_input.chr.par['amplitude_low_limit'])
        self.genetic_parameters.genes_amplitude_limits.append(user_input.chr.par['amplitude_hight_limit'])
        self.genetic_parameters.n_chromosomes = user_input.chr.par['n_chromosomes']



class Chromosome():
    def __init__(self):
        self.ampl_genes = None  # n_ampitudes
        self.J = None
        self.field = Field()
        self.prop_psi = prop.PropagatorEulero2Order()


    def init_chromosome(self, dt, molecule, field, pcm):
        self.field = field
        self.prop_psi.set_propagator(dt, molecule, pcm)
        self.ampl_genes = np.concatenate((self.field.parameters['fi'], self.field.parameters['fi_cos']), axis = 1)



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
        self.n_genes = None
        self.ampl_limits = []
        self.omegas_matrix = None


        self.chromosomes = []
        self.Js = []


    def iterate(self, current_iteration):
        if current_iteration != 0:
            self.evolve()
            for i in range(self.n_chromosomes):
                self.chromosomes[i].prop_psi.propagate_n_step(self.nstep, self.chromosomes[i].field.field)
            self.calc_J()
            self.check_convergence()
        else:
            for i in range(self.n_chromosomes):
                self.chromosomes[i].prop_psi.propagate_n_step(self.nstep, self.chromosomes[i].field.field)




    def init_oc_iterator(self, molecule, starting_field, pcm, oc_parameters, alpha_t):
        self.nstep = oc_parameters.nstep
        self.dt = oc_parameters.dt
        self.alpha_t = alpha_t
        self.target_state = oc_parameters.target_state
        #genetic
        self.n_genes = 2* starting_field.parameters['fi'].size()
        self.omegas_matrix = starting_field.parameters['omega']
        self.n_chromosomes = oc_parameters.iterator_parameters.n_chromosomes
        self.ampl_limits = oc_parameters.iterator_parameters.genes_amplitude_limits
        #genero i cromosomi con ampiezze random
        for i in range(self.n_chromosomes):
            chromosome = Chromosome()
            field = self.create_random_field(starting_field)
            chromosome.init_chromosome(self.dt, molecule, field, pcm)
            self.chromosomes.append(chromosome)


    def evolve(self):
        pass

    def init_evolutionary_algorithm(self):
        pass


    def create_random_field(self, starting_field):
        field = starting_field
        ampl = np.zeros([self.n_genes, 3])
        ampl[:, 0] = np.random.uniform((self.ampl_limits[0], self.ampl_limits[1], self.n_genes))
        ampl[:, 1] = np.random.uniform((self.ampl_limits[0], self.ampl_limits[1], self.n_genes))
        ampl[:, 2] = np.random.uniform((self.ampl_limits[0], self.ampl_limits[1], self.n_genes))
        field.parameters['fi'] = ampl[:self.n_genes/2]
        field.parameters['fi_cos'] = ampl[self.n_genes / 2:]
        #self.field.chose_field(self.field.field_type)
        field.chose_field('sum')
        return field



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
        #self.Js = []
        for i in range(self.n_chromosomes):
            #self.Js.append(
            self.chromosomes[i].J = np.real(af.projector_mean_value(self.chromosomes[i].prop_psi.mol.wf.ci, self.target_state) \
                         - af.alpha_field_J_integral(self.chromosomes[i].field.field, self.alpha_t, self.dt))

        self.chromosomes.sort(key=lambda x: x.J, reverse=True)
        #self.Js, self.chromosomes = (list(t) for t in zip(*sorted(zip(self.Js, self.chromosomes), reverse = True)))
        self.J = self.chromosomes[0].J



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
