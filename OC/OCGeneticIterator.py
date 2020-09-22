import time

import numpy as np
import pandas as pd
import multiprocessing as mp
import concurrent.futures


from parameters.GeneticParameters import GeneticParameters
from read_and_set.read import auxiliary_functions as af
from OC.ABCOCIterator import ABCOCIterator
from parameters.OCIteratorParameters import OCIteratorParameters
from SystemObj import DiscreteTimePar
from field.Field import Field, Func_tMatrix
from copy import deepcopy
from propagator import PropagatorsEulero as prop

from deap import base
from deap import creator
from deap import tools

import random

###################GENETIC ALGORITHMS#################
#sequential: first mate with crossover, then mutate all chromosomes

#mixed: mates n_tot - n_mutate chromosomes and mutate n_mutate chromosomes
#       only if mating and mutation are separated it is possible to evaluate mutation success

#restart: if stucked in local minima restart from a different point in the space

#clustering:

#dcrab:
#################SIGMA##########################
#adaptive_on_mutate_success:  depending on the muate successes, defined as eta > eta_thr, sigma changes
#                             only compatible with mixed genetic algorithm







Evolutionary_Algorithm_dict = {'cxUniform': tools.cxUniform,
                               'mutGaussian': tools.mutGaussian,
                               'selBest': tools.selBest}


class OCGeneticIterator(ABCOCIterator):
    def __init__(self):
        super().__init__()
        self.par = OCIteratorParameters()
        self.discrete_t_par = DiscreteTimePar()
        self.genetic_par = GeneticParameters()
        self.genetic_algorithms = base.Toolbox()

        self.field_psi_matrix = Func_tMatrix()
        self.psi_coeff_t_matrix = Func_tMatrix()
        self.dict_out = {}

        #genetic
        self.initial_c0 = None
        self.chromosomes = [] #can be a DEAP class or personal chromosome class
        self.n_mutation_succes = 0
        self.n_mutated = 0


    def iterate(self, current_iteration):
        print(current_iteration)
        self.evolve()
        self.check_convergence()

    def check_convergence(self):
        J_prev_tmp = np.copy(self.par.J)
        self.calc_J()
        self.par.convergence_t = self.par.J - J_prev_tmp
        self.par.convergence_t = self.par.convergence_t[0]

    def calc_J(self):
        self.chromosomes.sort(key=lambda x: x.J, reverse=True)
        self.par.J = self.chromosomes[0].J.values
        self.field_psi_matrix = self.chromosomes[0].field.field


    def init_output_dictionary(self):
        self.dict_out['log_file'] = self.get_log_file_out
        self.dict_out['final_pop'] = self.get_final_pop
        self.dict_out['pop_t'] = self.get_pop_t
        self.dict_out['field_t'] = self.get_field_t
        self.dict_out['field_ampl'] = self.get_field_ampl

    def init(self, molecule, starting_field, medium, alpha_t, oc_input, iterator_config_input):
        self.discrete_t_par.nstep = oc_input.nstep
        self.discrete_t_par.dt = oc_input.dt
        self.par.target_state = oc_input.target_state
        self.par.alpha_t = alpha_t
        self.par.convergence_t = 99999
        self.par.J = 99999

        self.field_psi_matrix = deepcopy(starting_field.field)

        #self.psi_coeff_t_matrix = np.zeros([self.simulation_par.nstep + 1, molecule.wf.n_ci], dtype=complex)
        self.init_output_dictionary()

        self.init_genetic(molecule, starting_field, medium, iterator_config_input)

    def init_genetic(self, molecule, starting_field, medium, genetic_input):
        self.genetic_par.genetic_algorithm = genetic_input.genetic_algorithm
        self.genetic_par.n_chromosomes = genetic_input.n_chromosomes
        self.genetic_par.n_selected_chr = genetic_input.n_selected_chr
        self.genetic_par.amplitude_lim = genetic_input.amplitude_lim

        self.genetic_par.mate = genetic_input.mate
        self.genetic_par.mate_probability = genetic_input.mate_probability

        self.genetic_par.mutate = genetic_input.mutate
        self.genetic_par.mutate_probability = genetic_input.mutate_probability
        self.genetic_par.n_mutate = genetic_input.n_mutate
        self.genetic_par.mutate_mu = genetic_input.mutate_mu
        self.genetic_par.mutate_starting_sigma = genetic_input.mutate_starting_sigma
        self.genetic_par.eta_thr = genetic_input.eta_thr
        self.genetic_par.q = genetic_input.q
        self.genetic_par.select = genetic_input.select
        self.genetic_par.n_amplitudes = starting_field.par.fi.size
        self.genetic_par.omegas_matrix = starting_field.par.omega


        self.init_chromosomes(molecule, starting_field, medium)
    #qui inizializzo l'algoritmo genetico che sto usando, che per adesso è solo uno
        self.init_evolutionary_algorithm(genetic_input)

    def init_chromosomes(self, molecule, starting_field, medium):
        self.init_DEAP_chromosomes(molecule, starting_field, medium)

    def create_random_ampl_rounded(self):
        number = round(random.uniform(-self.genetic_par.amplitude_lim, self.genetic_par.amplitude_lim),4)
        return number

    def init_DEAP_chromosomes(self, molecule, starting_field, medium):
        print("init_chromo")
        toolbox= base.Toolbox()
    # fitness è una classe astratta.
    # Creo creator.J, una classe figlia di fitness che ha pesi positivi quindi massimizza
        creator.create("J", base.Fitness, weights=(1.0,))
    # Definisco un oggetto che chiamo Chromosome, come attributi ha J, un oggetto campo e un oggetto propagator.
    # E poi di default ha un vettore da riempire che quando chiami il cromosoma senza attributi viene fuori quello
        creator.create("Chromosome", list, J=creator.J, field = Field(), prop_psi = prop.PropagatorEulero2Order())
    #creo un'istanza cromosoma singolo, di tipo Chromosome.
    # Creo n_amplitudes ampiezze con create_random_ampl_rounded e
    # automaticamente finiscono nel vettore vuoto di Chromosome
        toolbox.register("single_chromosome",
                         tools.initRepeat, creator.Chromosome,
                         self.create_random_ampl_rounded,
                         n=self.genetic_par.n_amplitudes)
    #faccio una n:cromosomi e li caccio in una lista
        toolbox.register('chromosomes_population',
                         tools.initRepeat,
                         list,
                         toolbox.single_chromosome)
        self.chromosomes = toolbox.chromosomes_population(n=self.genetic_par.n_chromosomes) #create all chromosomes
    #creo un oggetto prop_psi ....
        prop_psi = prop.PropagatorEulero2Order()
        prop_psi.set_propagator(molecule, medium)
        self.initial_c0 = prop_psi.mol.wf.ci
    #... e dentro a ogni cromosoma della mia lista inizializzo il campo, il propagatore e metto J = 0.5
    #Sono tutti uguali ma le ampiezze del mio vettore sono diverse perchè le ho fatte random prima
        for i in range(len(self.chromosomes)):
            self.chromosomes[i].field = deepcopy(starting_field)
            self.chromosomes[i].prop_psi = deepcopy(prop_psi)
            self.chromosomes[i].J.values = [0.5]
            if (starting_field.par.field_type == "restart_genetic"):
                self.assign_chromosome_values(self.chromosomes[i], starting_field.par.fi.reshape((1,-1)))
        self.par.J = [99999]


    def assign_chromosome_values(self, chromosome, reshaped_starting_field_fi):
        n = reshaped_starting_field_fi.size
        for j in range(n):
            chromosome[j] = reshaped_starting_field_fi[0][j] + random.uniform(-self.genetic_par.mutate_starting_sigma/10, -self.genetic_par.mutate_starting_sigma/10)


    def reset_chromosomes(self):
        for i in range(len(self.chromosomes)):
            self.chromosomes[i].prop_psi.medium.reset_medium(self.chromosomes[i].prop_psi.mol,
                                                             self.chromosomes[i].field.field)


    def evaluate_J_DEAP_chromosome(self, chro):
        chro.field.par.fi = np.asarray(chro).reshape((-1, 3))
        chro.field.chose_field('sum', discrete_t_par = self.discrete_t_par)
        chro.prop_psi.mol.wf.set_wf(self.initial_c0, 1)
        chro.prop_psi.propagate_n_step(self.discrete_t_par, chro.field.field)
        J = np.real(af.projector_mean_value(chro.prop_psi.mol.wf.ci,
                                            self.par.target_state)
                    - self.alpha_field_J_integral_chromosome(chro.field.field))
        success =  J - chro.J.values[0]
        if success != 0:
            self.n_mutated +=1
        if success > 0:
            self.n_mutation_succes +=1
        chro.J.values = [J]
        return chro


    def evolve(self):
        if self.genetic_par.genetic_algorithm == "mixed":
            self.evolve_mixed_DEAP()
        elif self.genetic_par.genetic_algorithm == "sequential":
            self.evolve_sequential_DEAP()

    def evolve_sequential_DEAP(self):
        n_children_each = self.genetic_par.n_chromosomes / self.genetic_par.n_selected_chr
        if n_children_each.is_integer() == False:
            n_children_each = int(n_children_each) + 1
        else:
            n_children_each = int(n_children_each)
        print("evolve")
        # Select the next generation individuals
        new = []
        selected = self.genetic_algorithms.select(self.chromosomes, fit_attr='J')
        selected = self.genetic_algorithms.clone(selected)
        # Apply crossover on the offspring
        for i in range(n_children_each):
            for first, second in zip(selected[::2], selected[1::2]):
                if random.random() < self.genetic_par.mate_probability:
                    a, b = self.genetic_algorithms.mate(deepcopy(first), deepcopy(second))
                else:
                    a, b = deepcopy(first), deepcopy(second)
                new.append(deepcopy(a))
                new.append(deepcopy(b))
            random.shuffle(selected)
        new = new[:self.genetic_par.n_chromosomes]
        # Apply mutation on the offspring
        for mutant in new:
            if random.random() < self.genetic_par.mutate_probability:
                self.genetic_algorithms.mutate(mutant)
        self.check_bounds_matrix(new)
        with concurrent.futures.ProcessPoolExecutor() as executor:
            new = list(executor.map(self.genetic_algorithms.evaluate, new))
        self.chromosomes[:] = new



    #quello del paper, con parte crossover e parte mutati
    def evolve_mixed_DEAP(self):
    #fa la stessa cosa che nel sequential, fa crossover con tutti. La valutazione invece la fa solo nel sottoset dei mutati
        n_mate_children_each = self.genetic_par.n_chromosomes / self.genetic_par.n_selected_chr
        if n_mate_children_each.is_integer() == False:
            n_mate_children_each = int(n_mate_children_each) + 1
        else:
            n_mate_children_each = int(n_mate_children_each)
        print("evolve")
        # Select the next generation individuals
        new_crossover = []
    #l'algoritmo è già definito. Quindi seleziono da cromosomi i migliori secondo J
        selected = self.genetic_algorithms.select(self.chromosomes, fit_attr='J')
        selected = self.genetic_algorithms.clone(selected)
    # Apply crossover on the offspring
    # li prendo tutti e ne genero quanti me ne servono mescolandoli
        for i in range(n_mate_children_each):
    #se tolgo questo il primo abbinamento me lo fa tra i vicini di J, se lo lascio è random
            random.shuffle(selected)
            for first, second in zip(selected[::2], selected[1::2]):
                if random.random() < self.genetic_par.mate_probability:
                    a, b = self.genetic_algorithms.mate(deepcopy(first), deepcopy(second))
                else:
                    a, b = deepcopy(first), deepcopy(second)
                new_crossover.append(deepcopy(a))
                new_crossover.append(deepcopy(b))

        len_cross = self.genetic_par.n_chromosomes - self.genetic_par.n_mutate
        new_crossover = new_crossover[:len_cross]

    # Apply mutation on the offspring
    # solo il sottoset

        new_mutation = []
        print("sigma:")
        print(self.genetic_algorithms.mutate.keywords['sigma'])
        for i in range(self.genetic_par.n_mutate):
            new_mutation.append(deepcopy(selected[i]))
            if random.random() < self.genetic_par.mutate_probability:
                self.genetic_algorithms.mutate(new_mutation[-1])
        self.check_bounds_matrix(new_mutation)
        #calcolo la J del crossover e la J della mutazione. in entrambi i casi valuta il numero di successi
        # ma resetto prima della mutazione perchè mi interessano solo quelli della mutazione
        new_crossover = list(self.genetic_algorithms.map(self.genetic_algorithms.evaluate, new_crossover))
        self.n_mutation_succes = 0
        self.n_mutated = 0
        new_mutation = list(self.genetic_algorithms.map(self.genetic_algorithms.evaluate, new_mutation))
        #updato sigma della mutazione secondo il successo di questa generazione
        self.update_sigma_mutation()
        #lo sostituisco alla mia popolazione
        new = new_crossover + new_mutation
        self.chromosomes[:] = new


    def update_sigma_mutation(self):
        print("mutation succes " + str(self.n_mutation_succes) + "n mutated " + str(self.n_mutated))
        eta = self.n_mutation_succes/(self.n_mutated)
        print("eta: " + str(eta))
        if eta < self.genetic_par.eta_thr:
            self.genetic_par.mutate_starting_sigma = self.genetic_par.mutate_starting_sigma*self.genetic_par.q
        else:
            self.genetic_par.mutate_starting_sigma = self.genetic_par.mutate_starting_sigma / self.genetic_par.q
        self.genetic_algorithms.mutate.keywords['sigma'] = self.genetic_par.mutate_starting_sigma
        self.n_mutation_succes = 0

    def check_bounds_matrix(self, matrix):
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                if matrix[i][j] > self.genetic_par.amplitude_lim:
                    matrix[i][j] = self.genetic_par.amplitude_lim
                elif matrix[i][j] < -self.genetic_par.amplitude_lim:
                    matrix[i][j] = -self.genetic_par.amplitude_lim


    def check_bounds_vector(self, vector):
        for i in range(len(vector)):
            if vector[i] > self.genetic_par.amplitude_lim:
                vector[i] = self.genetic_par.amplitude_lim
            elif vector[i] < -self.genetic_par.amplitude_lim:
                vector[i] = -self.genetic_par.amplitude_lim


    def init_evolutionary_algorithm(self, iterator_parameters):
        self.genetic_algorithms.register("mate", Evolutionary_Algorithm_dict[iterator_parameters.mate],
                                         indpb=self.genetic_par.mate_probability)
        self.genetic_algorithms.register("mutate",
                                         Evolutionary_Algorithm_dict[iterator_parameters.mutate],
                                         mu=self.genetic_par.mutate_mu,
                                         sigma=self.genetic_par.mutate_starting_sigma,
                                         indpb=self.genetic_par.mutate_probability)
        self.genetic_algorithms.register("select",
                                         Evolutionary_Algorithm_dict[iterator_parameters.select],
                                         k=self.genetic_par.n_selected_chr)
        self.reset_chromosomes()
        self.genetic_algorithms.register('evaluate',
                                         self.evaluate_J_DEAP_chromosome)

    def get_log_file_out(self):
        integral_field = np.real(self.field_J_integral())
        norm_proj = np.real(af.projector_mean_value(self.chromosomes[0].prop_psi.mol.wf.ci,
                                                    self.par.target_state) /
                            (np.dot(self.chromosomes[0].prop_psi.mol.wf.ci,
                                    np.conj(self.chromosomes[0].prop_psi.mol.wf.ci))))
        return np.array([self.par.convergence_t, self.par.J[0], norm_proj, integral_field])


    def get_final_pop(self):
        final_pop = np.real(af.population_from_wf_vector(self.chromosomes[0].prop_psi.mol.wf.ci))
        return final_pop


    def get_pop_t(self):
        self.chromosomes[0].prop_psi.mol.wf.set_wf(self.initial_c0, 1)
        psi_coeff_t_matrix = self.chromosomes[0].prop_psi.propagate_n_step(self.discrete_t_par,
                                                                           self.chromosomes[0].field.field)
        psi_coeff_t_matrix = np.insert(psi_coeff_t_matrix.f_xyz, 0, psi_coeff_t_matrix.time_axis, axis = 1)
        pop_t_matrix = np.real(af.population_from_wf_matrix(psi_coeff_t_matrix))
        return pop_t_matrix

    def get_field_t(self):
        field_t_matrix = np.insert(self.field_psi_matrix.f_xyz, 0, self.field_psi_matrix.time_axis, axis = 1)
        return field_t_matrix


    def get_field_ampl(self):
        return self.chromosomes[0].field.par.fi
        #return np.hstack((self.chromosomes[0].field.par.fi, self.chromosomes[0].field.par.fi_cos))

    def get_restart(self):
        separate1 = np.array([["###","omegas","###"],["###","###","###"]])
        separate2 = np.array([["###", "fi", "###"],["###","###","###"]])
        out = pd.DataFrame(separate1)
        out_tmp = pd.DataFrame(self.genetic_par.omegas_matrix)
        out = out.append(out_tmp)
        out_tmp = pd.DataFrame(separate2)
        out = out.append(out_tmp)
        out_tmp = pd.DataFrame(self.chromosomes[0].field.par.fi)
        out = out.append(out_tmp)
        return out

    def alpha_field_J_integral_chromosome(self, field):
        ax_square= field.f_xyz.ndim - 1
        ax_integral= field.f_xyz.ndim - 2
        f_square = np.sum(field.f_xyz * field.f_xyz, axis=ax_square) * self.par.alpha_t
        f_integral = np.sum(f_square, axis=ax_integral)
        out_integral = f_integral*self.discrete_t_par.dt
        return out_integral

