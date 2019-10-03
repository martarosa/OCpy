import numpy as np
import array

from read import auxiliary_functions as af
from OCIterator import OCIterator
from Chromosome import Chromosome
from field.Field import Field
from copy import deepcopy
import DeapWrapper as deap
from propagator import PropagatorsEulero as prop

from deap import base
from deap import creator
from deap import tools

import random

###################DEAP INITIALIZATION#################
#mate
#cxUniform
prob_mate = 0.5
#mutate
#mutGaussian
mu = 0.1
sigma = 0.1
prob_mutate = 0.1
#######################################################


DEAP_EA_dict = {'cxUniform': tools.cxUniform,
                'mutGaussian': tools.mutGaussian,
                'selbest': tools.selBest}




class OCDeapIterator():
    def __init__(self):
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

        #genetic
        self.n_amplitudes = None
        self.n_chromosomes = None
        self.n_evolved_chr = None
        self.ampl_min = None
        self.ampl_max = None
        self.omegas_matrix = None
        #DEAP
        self.toolbox = base.Toolbox()
        self.chromosomes =  None



        def iterate(self, current_iteration):
            pass

        def check_convergence(self):
            pass

        def calc_J(self):
            pass

        def init_output_dictionary(self):
            pass



    def init(self, oc_parameters, iterator_parameters, molecule, starting_field, pcm, alpha_t):
        self.nstep = oc_parameters.nstep
        self.dt = oc_parameters.dt
        self.target_state = oc_parameters.target_state
        self.alpha_t = alpha_t

        self.n_chromosomes = iterator_parameters.n_chromosomes
        self.n_amplitudes = 2 * starting_field.parameters['fi'].shape[0]
        self.omegas_matrix = starting_field.parameters['omega']
        self.ampl_min = iterator_parameters.amplitude_min
        self.ampl_max = iterator_parameters.amplitude_max
        self.n_evolved_chr = iterator_parameters.n_evolved_chr

        creator.create("J", base.Fitness, weights=(1.0,))
        creator.create("Chromosome", list, fitness=creator.J, field = Field(), propagator = prop.PropagatorEulero2Order())

        self.toolbox.register("create_random_ampl", random.uniform, self.ampl_min, self.ampl_max)
        self.toolbox.register("single_chromosome", tools.initRepeat, creator.Chromosome, self.toolbox.create_random_ampl, n=self.n_amplitudes)
        self.toolbox.register('chromosomes_population', tools.initRepeat, list, self.toolbox.single_chromosome)

        self.chromosomes = self.toolbox.chromosomes_population(n=self.n_chromosomes)

        prop_psi = prop.PropagatorEulero2Order()
        prop_psi.set_propagator(self.dt, molecule, pcm)
        for i in range(len(self.chromosomes)):
            self.chromosomes[i].field = deepcopy(starting_field)
            self.chromosomes[i].propagator = deepcopy(prop_psi)


    def amplitudes_to_field(self, i):
        self.chromosomes[i].field.parameters['fi'] = np.asarray(self.chromosomes[i][:int(self.n_amplitudes / 2)]).reshape(
            (-1, 3))
        self.chromosomes[i].field.parameters['fi_cos'] = np.asarray(self.chromosomes[i][:int(self.n_amplitudes / 2)]).reshape(
            (-1, 3))
        self.chromosomes[i].field.chose_field('sum')

    def evaluate_chromosomes(self, nstep):
        out_propagation = []
        for i in range(len(self.chromosomes)):
            self.amplitudes_to_field(i)
            out_propagation.append(self.chromosomes[i].propagator.propagate_n_step(nstep, self.chromosomes[i].field))
            self.chromosomes[i].fitness.values = self.calc_J_cromosome(i)

    def calc_J_chromosome(self, i):
        J = np.real(af.projector_mean_value(self.chromosomes[i].prop_psi.propagator_terms.mol.wf.ci,
                                            self.target_state)
                    - af.alpha_field_J_integral(self.chromosomes[i].field.field, self.alpha_t, self.dt))
        return J

    def init_DEAP_evolutionary_algorithm(self, iterator_parameters):
        if(iterator_parameters.mate == 'cxUniform'):
            indpb = [prob_mate] * self.n_amplitudes
            self.toolbox.register("mate", DEAP_EA_dict[iterator_parameters.mate], indpb = indpb)
        if(iterator_parameters.mutate == 'mutGaussian'):
            self.toolbox.register("mutate", DEAP_EA_dict[iterator_parameters.mutate], mu = mu, sigma = sigma, indpb = prob_mutate)
        if(iterator_parameters.select == 'selBest')
            self.toolbox.register("select", DEAP_EA_dict[iterator_parameters.select], self.n_evolved_chr)

        self.toolbox.register('evaluate', self.evaluate_chromosomes, nstep = self.nstep)







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



        #genetic
        self.n_amplitudes = None
        self.n_chromosomes = None
        self.n_evolved_chr = None
        self.ampl_min = None
        self.ampl_max = None
        self.omegas_matrix = None

        self.chromosomes =  []



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



    def evolve(self):
        #evolvere le ampiezze nei cromosomi e passare da ampiezze a campi
        pass




    def init_evolutionary_algorithm(self):
        pass

    def init(self, oc_parameters, iterator_parameters, molecule, starting_field, pcm, alpha_t):
        self.nstep = oc_parameters.nstep
        self.dt = oc_parameters.dt
        self.alpha_t = alpha_t
        self.target_state = oc_parameters.target_state
        #iterator_parameters = genetic
        self.n_amplitudes = 2 * starting_field.parameters['fi'].shape[0]
        self.omegas_matrix = starting_field.parameters['omega']
        self.n_chromosomes = iterator_parameters.n_chromosomes
        self.ampl_min = iterator_parameters.amplitude_min
        self.ampl_max = iterator_parameters.amplitude_max
        self.n_evolved_chr = iterator_parameters.n_evolved_chr
        #genero i cromosomi con ampiezze random
        for i in range(self.n_chromosomes):
            chromosome = Chromosome()
            field = self.create_random_starting_chromosome_field(starting_field)
            chromosome.init_chromosome(self.dt, molecule, field, pcm)
            self.chromosomes.append(deepcopy(chromosome))


    def create_random_starting_chromosome_field(self, starting_field):
        field = starting_field
        ampl = np.zeros([self.n_amplitudes, 3])
        ampl[:, 0] = np.random.uniform(self.ampl_min, self.ampl_max, self.n_amplitudes)
        ampl[:, 1] = np.random.uniform(self.ampl_min, self.ampl_max, self.n_amplitudes)
        ampl[:, 2] = np.random.uniform(self.ampl_min, self.ampl_max, self.n_amplitudes)
        field.parameters['fi'] = ampl[:int(self.n_amplitudes/2)]
        field.parameters['fi_cos'] = ampl[int(self.n_amplitudes / 2):]
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
        for i in range(self.n_chromosomes):
            #self.Js.append(
            self.chromosomes[i].J = np.real(af.projector_mean_value(
                self.chromosomes[i].prop_psi.propagator_terms.mol.wf.ci,
                self.target_state)
                                            - af.alpha_field_J_integral(
                self.chromosomes[i].field.field,
                self.alpha_t,
                self.dt))

        self.chromosomes.sort(key=lambda x: x.J, reverse=True)
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





