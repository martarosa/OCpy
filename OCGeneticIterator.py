import numpy as np
import array
import sys

from read import auxiliary_functions as af
from OCIterator import OCIterator, OCIteratorAttr
from Chromosome import Chromosome
from field.Field import Field
from copy import deepcopy
from propagator import PropagatorsEulero as prop

from deap import base
from deap import creator
from deap import tools

import random

###################DEAP INITIALIZATION#################
#mate
mate_probability = 1 #given 2 list which is the probability they mate
#cxUniform
cx_uniform_probability = 0.5
#mutate
mutate_probability = 0.2
n_mutate = 2
#mutGaussian
mu = 0
eta_thr = 0.1
q = 0.2
#######################################################


Evolutionary_Algorithm_dict = {'DEAP_cxUniform': tools.cxUniform,
                               'DEAP_mutGaussian': tools.mutGaussian,
                               'DEAP_selBest': tools.selBest}



class OCGeneticIterator(OCIterator):
    def __init__(self):
        super().__init__()
        self.class_attributes = OCIteratorAttr()

        self.chromosomes = [] #can be a DEAP class or personal chromosome class

        #genetic
        self.initial_c0 = None

        self.n_amplitudes = None
        self.n_chromosomes = None
        self.n_selected_chr = None
        self.ampl_min = None
        self.ampl_max = None
        self.omegas_matrix = None

        self.n_mutation_succes = None
        self.mutate_sigma = 0.01

        self.deap = None
        self.evolutionary_algorithms = None


    def iterate(self, current_iteration):
        print(current_iteration)
        self.evolve()
        self.check_convergence()


    def init_output_dictionary(self):
        self.class_attributes.dict_out['log_file'] = self.get_log_file_out
        self.class_attributes.dict_out['final_pop'] = self.get_final_pop
        self.class_attributes.dict_out['pop_t'] = self.get_pop_t
        self.class_attributes.dict_out['field_t'] = self.get_field_t
        self.class_attributes.dict_out['field_ampl'] = self.get_field_ampl


    def init(self, oc_parameters, iterator_parameters, molecule, starting_field, pcm, alpha_t):
        self.class_attributes.nstep = oc_parameters.nstep
        self.class_attributes.dt = oc_parameters.dt
        self.class_attributes.alpha_t = alpha_t
        self.class_attributes.target_state = oc_parameters.target_state
        self.class_attributes.convergence_t = 99999
        self.class_attributes.J = 99999
        self.class_attributes.field_psi_matrix = np.zeros([self.class_attributes.nstep, 3])
        self.class_attributes.psi_coeff_t = np.zeros([self.class_attributes.nstep + 1, molecule.wf.n_ci], dtype=complex)
        self.init_output_dictionary()

        #iterator_parameters = genetic
        self.n_chromosomes = iterator_parameters.n_chromosomes
        self.n_amplitudes = starting_field.parameters['fi'].size
        self.omegas_matrix = starting_field.parameters['omega']
        self.ampl_min = iterator_parameters.amplitude_min
        self.ampl_max = iterator_parameters.amplitude_max
        self.n_selected_chr = iterator_parameters.n_evolved_chr
        self.deap = iterator_parameters.DEAP

        self.init_chromosomes(molecule, starting_field, pcm)
        self.init_evolutionary_algorithm(iterator_parameters)


    def calc_J(self):
        self.chromosomes.sort(key=lambda x: x.J, reverse=True)
        self.class_attributes.J = self.chromosomes[0].J.values
        self.class_attributes.field_psi_matrix = self.chromosomes[0].field.field


    def check_convergence(self):
        J_prev_tmp = np.copy(self.class_attributes.J)
        self.calc_J()
        self.class_attributes.convergence_t = self.class_attributes.J - J_prev_tmp
        self.class_attributes.convergence_t = self.class_attributes.convergence_t[0]



    def init_chromosomes(self, molecule, starting_field, pcm):
        if self.deap == 'true':
            self.init_DEAP_chromosomes(molecule, starting_field, pcm)
        else:
            self.init_general_chromosome(molecule, starting_field, pcm)


    def create_random_ampl_rounded(self):
        number = round(random.uniform(self.ampl_min, self.ampl_max),4)
        return number

    def init_DEAP_chromosomes(self, molecule, starting_field, pcm):
        print("init_chromo")
        toolbox= base.Toolbox()
        creator.create("J", base.Fitness, weights=(1.0,))
        creator.create("Chromosome", list, J=creator.J, field = Field(), prop_psi = prop.PropagatorEulero2Order())
        #toolbox.register("create_random_ampl", random.uniform, self.ampl_min, self.ampl_max)
        #toolbox.register("create_random_ampl_rounded", round, toolbox.create_random_ampl, ndigits=4)
        toolbox.register("single_chromosome", tools.initRepeat, creator.Chromosome, self.create_random_ampl_rounded, n=self.n_amplitudes)
        toolbox.register('chromosomes_population', tools.initRepeat, list, toolbox.single_chromosome)
        self.chromosomes = toolbox.chromosomes_population(n=self.n_chromosomes) #create all chromosomes
        prop_psi = prop.PropagatorEulero2Order()
        prop_psi.set_propagator(self.class_attributes.dt, molecule, pcm)
        self.initial_c0 = prop_psi.propagator_terms.mol.wf.ci
        for i in range(len(self.chromosomes)):
            self.chromosomes[i].field = deepcopy(starting_field)
            self.chromosomes[i].prop_psi = deepcopy(prop_psi)
            self.chromosomes[i].J.values = [0.5]


    def init_general_chromosome(self, molecule, starting_field, pcm):
        chromosome = Chromosome()
        chromosome.init_chromosome(self.class_attributes.dt, molecule, starting_field, pcm)
        self.chromosomes = self.n_chromosomes * chromosome
        for i in range(self.n_chromosomes):
            rand = []
            for j in range(self.n_amplitudes):
                rand.append(random.uniform(self.ampl_min, self.ampl_max))
            self.chromosomes[i].amplitudes = deepcopy(rand)
            self.chromosomes[i].J = [0.5]


    def evaluate_J_general_chromosomes(self):
        for i in range(len(self.chromosomes)):
            self.chromosomes[i].amplitudes_to_field()
            self.chromosomes[i].prop_psi.propagate_n_step(self.class_attributes.nstep, self.chromosomes[i].field.field)
            self.chromosomes[i].calc_J(self.class_attributes.target_state,
                                       self.class_attributes.alpha_t,
                                       self.class_attributes.dt)


    def evaluate_J_DEAP_chromosome(self, chro):
        chro.field.parameters['fi'] = np.asarray(chro).reshape((-1, 3))
        chro.field.chose_field('sum')
        #chro.field.chose_field('sum_pip')
        chro.prop_psi.propagator_terms.mol.wf.set_wf(self.initial_c0, 1)
        chro.prop_psi.propagate_n_step(self.class_attributes.nstep, chro.field.field)
        pop= np.real(af.projector_mean_value(chro.prop_psi.propagator_terms.mol.wf.ci,
                                            self.class_attributes.target_state))
        field = np.real(af.alpha_field_J_integral(chro.field.field,
                                                self.class_attributes.alpha_t,
                                                self.class_attributes.dt))

        J = np.real(af.projector_mean_value(chro.prop_psi.propagator_terms.mol.wf.ci,
                                            self.class_attributes.target_state)
                    - af.alpha_field_J_integral(chro.field.field,
                                                self.class_attributes.alpha_t,
                                                self.class_attributes.dt))
        if J > chro.J.values:
            self.n_mutation_succes += 1
        return [J]


    def evolve(self):
        if self.deap == 'true':
            self.evolve_mixed_DEAP()


    def evolve_general(self):
        # Select the next generation individuals
        selected = self.evolutionary_algorithms.select(self.chromosomes, fit_attr='J')

        new = []
        # Clone the selected individuals
        selected = self.evolutionary_algorithms.clone(selected)
        # Apply crossover on the offspring
        for i in range(int(self.n_chromosomes/self.n_selected_chr) + 1):
            for first, second in zip(selected[::2], selected[1::2]):
                if random.random() < mate_probability:
                    a, b = self.evolutionary_algorithms.mate(deepcopy(first.amplitudes), deepcopy(second.amplitudes))
                    new.append(deepcopy(a))
                    new.append(deepcopy(b))
            random.shuffle(selected)
        new = new[:self.n_chromosomes]

        # Apply mutation on the offspring
        for mutant in new:
            if random.random() < mutate_probability:
                self.evolutionary_algorithms.mutate(mutant)
        self.check_bounds(new)
        for i in range(self.n_chromosomes):
            self.chromosomes[i].amplitudes = new[i]
        self.evaluate_J_general_chromosomes()





    def evolve_subsequent_DEAP(self):
        n_children_each = self.n_chromosomes / self.n_selected_chr
        if n_children_each.is_integer() == False:
            n_children_each = int(n_children_each) + 1
        else:
            n_children_each = int(n_children_each)
        print("evolve")
        # Select the next generation individuals
        new = []
        selected = self.evolutionary_algorithms.select(self.chromosomes, fit_attr='J')
        selected = self.evolutionary_algorithms.clone(selected)
        # Apply crossover on the offspring
        for i in range(n_children_each):
            for first, second in zip(selected[::2], selected[1::2]):
                if random.random() < mate_probability:
                    a, b = self.evolutionary_algorithms.mate(deepcopy(first), deepcopy(second))
                else:
                    a, b = deepcopy(first), deepcopy(second)
                new.append(deepcopy(a))
                new.append(deepcopy(b))

            random.shuffle(selected)
        new = new[:self.n_chromosomes]

        # Apply mutation on the offspring
        for mutant in new:
            if random.random() < mutate_probability:
                self.evolutionary_algorithms.mutate(mutant)
        self.check_bounds(new)
        J = self.evolutionary_algorithms.map(self.evolutionary_algorithms.evaluate, new)
        self.update_sigma_mutation(self.n_chromosomes)
        for ind, fit in zip(new, J):
            ind.J.values = fit
        self.chromosomes[:] = new



    def evolve_mixed_DEAP(self):
        n_children_each = self.n_chromosomes/self.n_selected_chr

        if n_children_each.is_integer() == False:
            sys.exit("wrong numer of chromosome to evolve, n_chromosomes/n_selected must be integer")
        print("evolve")
        # Select the next generation individuals
        new_crossover = []
        selected = self.evolutionary_algorithms.select(self.chromosomes, fit_attr='J')
        selected = self.evolutionary_algorithms.clone(selected)
        # Apply crossover on the offspring
        for i in range(int(n_children_each-n_mutate)):
            for first, second in zip(selected[::2], selected[1::2]):
                if random.random() < mate_probability:
                    a, b = self.evolutionary_algorithms.mate(deepcopy(first), deepcopy(second))
                else:
                    a, b = deepcopy(first), deepcopy(second)
                new_crossover.append(deepcopy(a))
                new_crossover.append(deepcopy(b))
            random.shuffle(selected)

        # Apply mutation on the offspring
        new_mutation = []
        for i in range(n_mutate):
            for j in range(self.n_selected_chr):
                new_mutation.append(deepcopy(selected[j]))
                if random.random() < mutate_probability:
                    self.evolutionary_algorithms.mutate(new_mutation[-1])
        self.check_bounds(new_mutation)

        J_crossover = list(self.evolutionary_algorithms.map(self.evolutionary_algorithms.evaluate, new_crossover))
        self.n_mutation_succes = 0
        J_mutation = list(self.evolutionary_algorithms.map(self.evolutionary_algorithms.evaluate, new_mutation))
        self.update_sigma_mutation(n_mutate * self.n_selected_chr)
        new = new_crossover + new_mutation
        J = J_crossover + J_mutation
        for ind, fit in zip(new, J):
            ind.J.values = fit
        self.chromosomes[:] = new


    def update_sigma_mutation(self, n_mutated):
        eta = self.n_mutation_succes/n_mutated
        if eta <= eta_thr:
            self.mutate_sigma = self.mutate_sigma*q
        else:
            self.mutate_sigma = self.mutate_sigma / q

    def check_bounds(self, matrix):
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                if matrix[i][j] > self.ampl_max:
                    matrix[i][j] = self.ampl_max
                elif matrix[i][j] < self.ampl_min:
                    matrix[i][j] = self.ampl_min




    def init_evolutionary_algorithm(self, iterator_patameters):
        self.evolutionary_algorithms = base.Toolbox()
        self.evolutionary_algorithms.register("mate", Evolutionary_Algorithm_dict[iterator_patameters.mate], indpb=cx_uniform_probability)
        self.evolutionary_algorithms.register("mutate", Evolutionary_Algorithm_dict[iterator_patameters.mutate], mu=mu, sigma=self.mutate_sigma, indpb=mutate_probability)
        self.evolutionary_algorithms.register("select", Evolutionary_Algorithm_dict[iterator_patameters.select], k=self.n_selected_chr)
        self.evolutionary_algorithms.register('evaluate', self.evaluate_J_DEAP_chromosome)




    def get_log_file_out(self):
        integral_field = np.real(af.field_J_integral(self.class_attributes.field_psi_matrix, self.class_attributes.dt))
        norm_proj = np.real(af.projector_mean_value(self.chromosomes[0].prop_psi.propagator_terms.mol.wf.ci,
                                                    self.class_attributes.target_state)/
                            (np.dot(self.chromosomes[0].prop_psi.propagator_terms.mol.wf.ci,
                                    np.conj(self.chromosomes[0].prop_psi.propagator_terms.mol.wf.ci))))
        return np.array([self.class_attributes.convergence_t, self.class_attributes.J[0], norm_proj, integral_field])


    def get_final_pop(self):
        final_pop = np.real(af.population_from_wf_vector(self.chromosomes[0].prop_psi.propagator_terms.mol.wf.ci))
        return final_pop


    def get_pop_t(self):
        self.chromosomes[0].prop_psi.propagator_terms.mol.wf.set_wf(self.initial_c0, 1)
        self.class_attributes.psi_coeff_t = self.chromosomes[0].prop_psi.propagate_n_step(self.class_attributes.nstep,
                                                                                          self.chromosomes[0].field.field)
        pop_t = np.real(af.population_from_wf_matrix(self.class_attributes.psi_coeff_t))
        return pop_t


    def get_field_t(self):
        return self.class_attributes.field_psi_matrix


    def get_field_ampl(self):
        return self.chromosomes[0].field.parameters['fi']

    def get_restart(self):
        out = np.concatenate((self.omegas_matrix, self.chromosomes[0].field.parameters['fi']))
        return out











