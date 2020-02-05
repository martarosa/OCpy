import numpy as np
import array
import sys

from read import auxiliary_functions as af
from ABCOCIterator import ABCOCIterator, OCIteratorParameters
from Chromosome import Chromosome
from field.Field import Field
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



#mate
mate_probability = 1 #given 2 list which is the probability they mate
#cxUniform
cx_uniform_probability = 0.5
#mutate
mutate_probability = 0.2
n_mutate = 6
#mutGaussian
mu = 0.00
eta_thr = 0.6
q = 0.9

#######################################################




Evolutionary_Algorithm_dict = {'DEAP_cxUniform': tools.cxUniform,
                               'DEAP_mutGaussian': tools.mutGaussian,
                               'DEAP_selBest': tools.selBest}



class OCGeneticIterator(ABCOCIterator):
    def __init__(self):
        super().__init__()
        self.oc_iterator_parameters = OCIteratorParameters()

        #genetic
        self.initial_c0 = None

        self.chromosomes = [] #can be a DEAP class or personal chromosome class
        self.n_chromosomes = None
        self.n_selected_chr = None

        self.n_amplitudes = None
        self.ampl_min = None
        self.ampl_max = None

        self.omegas_matrix = None

        self.n_mutation_succes = 0
        self.mutate_sigma = 0.01


        self.n_mutate = n_mutate

        self.deap = None
        self.evolutionary_algorithms = None


    def iterate(self, current_iteration):
        print(current_iteration)
        self.evolve()
        self.check_convergence()


    def init_output_dictionary(self):
        self.oc_iterator_parameters.dict_out['log_file'] = self.get_log_file_out
        self.oc_iterator_parameters.dict_out['final_pop'] = self.get_final_pop
        self.oc_iterator_parameters.dict_out['pop_t'] = self.get_pop_t
        self.oc_iterator_parameters.dict_out['field_t'] = self.get_field_t
        self.oc_iterator_parameters.dict_out['field_ampl'] = self.get_field_ampl


    def init(self, oc_input, molecule, starting_field, pcm, alpha_t):
        self.oc_iterator_parameters.nstep = oc_input.nstep
        self.oc_iterator_parameters.dt = oc_input.dt
        self.oc_iterator_parameters.alpha_t = alpha_t
        self.oc_iterator_parameters.target_state = oc_input.target_state
        self.oc_iterator_parameters.convergence_t = 99999
        self.oc_iterator_parameters.J = 99999
        self.oc_iterator_parameters.field_psi_matrix = np.zeros([self.oc_iterator_parameters.nstep, 3])
        self.oc_iterator_parameters.psi_coeff_t = np.zeros([self.oc_iterator_parameters.nstep + 1, molecule.wf.n_ci], dtype=complex)
        self.init_output_dictionary()

        #iterator_parameters = genetic
        self.n_chromosomes = oc_input.n_chromosomes
        self.n_amplitudes = starting_field.parameters['fi'].size
        #self.n_chro_generated_mutation = n_mutate * self.n_chromosomes
        self.omegas_matrix = starting_field.parameters['omega']
        self.ampl_min = oc_input.amplitude_min
        self.ampl_max = oc_input.amplitude_max
        self.n_selected_chr = oc_input.n_evolved_chr
        self.deap = oc_input.DEAP

        self.init_chromosomes(molecule, starting_field, pcm)
    #qui inizializzo l'algoritmo genetico che sto usando, che per adesso è solo uno
        self.init_evolutionary_algorithm(oc_input)


    def calc_J(self):
        self.chromosomes.sort(key=lambda x: x.J, reverse=True)
        self.oc_iterator_parameters.J = self.chromosomes[0].J.values
        self.oc_iterator_parameters.field_psi_matrix = self.chromosomes[0].field.field


    def check_convergence(self):
        J_prev_tmp = np.copy(self.oc_iterator_parameters.J)
        self.calc_J()
        self.oc_iterator_parameters.convergence_t = self.oc_iterator_parameters.J - J_prev_tmp
        self.oc_iterator_parameters.convergence_t = self.oc_iterator_parameters.convergence_t[0]


    def init_chromosomes(self, molecule, starting_field, pcm):
        if self.deap == 'true':
            self.init_DEAP_chromosomes(molecule, starting_field, pcm)
        else:
            self.init_general_chromosome(molecule, starting_field, pcm)


    def create_random_ampl_rounded(self):
        number = round(random.uniform(-0.005, 0.005),4)
        return number



    def init_DEAP_chromosomes(self, molecule, starting_field, pcm):
        print("init_chromo")
        toolbox= base.Toolbox()
    # fitness è una classe astratta.
    # Creo creator.J, una classe figlia di fitness che ha pesi positivi quindi massimizza
        creator.create("J", base.Fitness, weights=(1.0,))
    # Definisco un oggetto che chiamo Chromosome, come attributi ha J, un oggetto campo e un oggetto propagator.
    # E poi di default ha un vettore da riempire che quando chiami il cromosoma senza attributi viene fuori quello
    #def Chromosome()
    #    self. = []
    #    self.field = Field()
    #    self.prop_psi = prop.PropagatorEulero2Order()
        creator.create("Chromosome", list, J=creator.J, field = Field(), prop_psi = prop.PropagatorEulero2Order())
        #toolbox.register("create_random_ampl", random.uniform, self.ampl_min, self.ampl_max)
        #toolbox.register("create_random_ampl_rounded", round, toolbox.create_random_ampl, ndigits=4)

    #creo un'istanza cromosoma singolo, di tipo Chromosome.
    # Creo n_amplitudes ampiezze con create_random_ampl_rounded e
    # automaticamente finiscono nel vettore vuoto di Chromosome
        toolbox.register("single_chromosome", tools.initRepeat, creator.Chromosome, self.create_random_ampl_rounded, n=self.n_amplitudes)
    #faccio una n:cromosomi e li caccio in una lista
        toolbox.register('chromosomes_population', tools.initRepeat, list, toolbox.single_chromosome)
        self.chromosomes = toolbox.chromosomes_population(n=self.n_chromosomes) #create all chromosomes

    #creo un oggetto prop_psi ....
        prop_psi = prop.PropagatorEulero2Order()
        prop_psi.set_propagator(self.oc_iterator_parameters.dt, molecule, pcm)
        self.initial_c0 = prop_psi.propagator_terms.mol.wf.ci
    #... e dentro a ogni cromosoma della mia lista inizializzo il campo, il propagatore e metto J = 0.5
    #Sono tutti uguali ma le ampiezze del mio vettore sono diverse perchè le ho fatte random prima
        for i in range(len(self.chromosomes)):
            self.chromosomes[i].field = deepcopy(starting_field)
            self.chromosomes[i].prop_psi = deepcopy(prop_psi)
            self.chromosomes[i].J.values = [0.5]
        self.oc_iterator_parameters.J = [99999]



    def evaluate_J_DEAP_chromosome(self, chro):
    # prende l'array di Chromosme e lo mette nel campo come ampiezze
        chro.field.parameters['fi'] = np.asarray(chro).reshape((-1, 3))
    #genera il campo con forma sum
        chro.field.chose_field('sum_pip')
        #chro.field.chose_field('sum_pip')
    #setto il propagatore e lo propago con il campo appena generato
        chro.prop_psi.propagator_terms.mol.wf.set_wf(self.initial_c0, 1)
        chro.prop_psi.propagate_n_step(self.oc_iterator_parameters.nstep, chro.field.field)
    #calcolo popolazione e integrale del campo, per debug
        pop= np.real(af.projector_mean_value(chro.prop_psi.propagator_terms.mol.wf.ci,
                                             self.oc_iterator_parameters.target_state))
        field = np.real(af.alpha_field_J_integral(chro.field.field,
                                                  self.oc_iterator_parameters.alpha_t,
                                                  self.oc_iterator_parameters.dt))
    #calcolo J
        J = np.real(af.projector_mean_value(chro.prop_psi.propagator_terms.mol.wf.ci,
                                            self.oc_iterator_parameters.target_state)
                    - af.alpha_field_J_integral(chro.field.field,
                                                self.oc_iterator_parameters.alpha_t,
                                                self.oc_iterator_parameters.dt))

        success =  J - chro.J.values[0]
        if success > 0:
            self.n_mutation_succes +=1
        return [J]


    def evolve(self):
        if self.deap == 'true':
            self.evolve_subsequent_DEAP()



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
        for ind, fit in zip(new, J):
            ind.J.values = fit
        self.chromosomes[:] = new


    #quello del paper, con parte crossover e parte mutati
    def evolve_mixed_DEAP(self):
        n_mate = self.n_chromosomes-(self.n_selected_chr*n_mutate)
        print("evolve")
        # Select the next generation individuals
        new_crossover = []
    #l'algoritmo è già definito. Quindi seleziono da cromosomi i migliori secondo J
        selected = self.evolutionary_algorithms.select(self.chromosomes, fit_attr='J')
        selected = self.evolutionary_algorithms.clone(selected)
    # Apply crossover on the offspring
    # li prendo tutti e ne genero quanti me ne servono mescolandoli
        for i in range(n_mate):
    #se tolgo questo il primo abbinamento me lo fa tra i vicini di J, se lo lascio è random
            random.shuffle(selected)
            for first, second in zip(selected[::2], selected[1::2]):
                if random.random() < mate_probability:
                    a, b = self.evolutionary_algorithms.mate(deepcopy(first), deepcopy(second))
                else:
                    a, b = deepcopy(first), deepcopy(second)
                new_crossover.append(deepcopy(a))
                new_crossover.append(deepcopy(b))
    # Apply mutation on the offspring
    # li prendo tutti in ordine e li muto le volte che mi serve
        new_mutation = []
        print("sigma:")
        print(self.evolutionary_algorithms.mutate.keywords['sigma'])
        for i in range(self.n_selected_chr):
            for j in range(self.n_mutate):
                new_mutation.append(deepcopy(selected[i]))
                self.evolutionary_algorithms.mutate(new_mutation[-1])
        #self.check_bounds(new_mutation)

        #calcolo la J del crossover e la J della mutazione
        J_crossover = list(self.evolutionary_algorithms.map(self.evolutionary_algorithms.evaluate, new_crossover))
        self.n_mutation_succes = 0
        J_mutation = list(self.evolutionary_algorithms.map(self.evolutionary_algorithms.evaluate, new_mutation))
        print("mutation_success:")
        print(self.n_mutation_succes)
        #updato sigma della mutazione secondo il successo di questa generazione

        self.update_sigma_mutation()
        #genero il nuovo pool, lo zippo con le sue J
        new = new_crossover + new_mutation
        J = J_crossover + J_mutation
        for ind, fit in zip(new, J):
            ind.J.values = fit
    #lo sostituisco all amia popolazione
        self.chromosomes[:] = new




    def update_sigma_mutation(self):
        print("eta:")
        eta = self.n_mutation_succes/(self.n_mutate*self.n_selected_chr)
        print(eta)
        #print(eta)
        #print(self.mutate_sigma)
        if eta <= eta_thr:
            self.mutate_sigma = self.mutate_sigma*q
        else:
            self.mutate_sigma = self.mutate_sigma / q
        self.evolutionary_algorithms.mutate.keywords['sigma'] = self.mutate_sigma
        self.n_mutation_succes = 0

    def check_bounds(self, matrix):
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                if matrix[i][j] > self.ampl_max:
                    matrix[i][j] = self.ampl_max
                elif matrix[i][j] < self.ampl_min:
                    matrix[i][j] = self.ampl_min
                #elif (matrix[i][j] < 0.0005 and matrix[i][j] > -0.0005):
                #    matrix[i][j] = 0.0




    def init_evolutionary_algorithm(self, iterator_parameters):
        self.evolutionary_algorithms = base.Toolbox()
        self.evolutionary_algorithms.register("mate", Evolutionary_Algorithm_dict[iterator_parameters.mate], indpb=cx_uniform_probability)
        self.evolutionary_algorithms.register("mutate", Evolutionary_Algorithm_dict[iterator_parameters.mutate], mu=mu, sigma=self.mutate_sigma, indpb=mutate_probability)
        self.evolutionary_algorithms.register("select", Evolutionary_Algorithm_dict[iterator_parameters.select], k=self.n_selected_chr)
        self.evolutionary_algorithms.register('evaluate', self.evaluate_J_DEAP_chromosome)




    def get_log_file_out(self):
        integral_field = np.real(af.field_J_integral(self.oc_iterator_parameters.field_psi_matrix, self.oc_iterator_parameters.dt))
        norm_proj = np.real(af.projector_mean_value(self.chromosomes[0].prop_psi.propagator_terms.mol.wf.ci,
                                                    self.oc_iterator_parameters.target_state) /
                            (np.dot(self.chromosomes[0].prop_psi.propagator_terms.mol.wf.ci,
                                    np.conj(self.chromosomes[0].prop_psi.propagator_terms.mol.wf.ci))))
        return np.array([self.oc_iterator_parameters.convergence_t, self.oc_iterator_parameters.J[0], norm_proj, integral_field])


    def get_final_pop(self):
        final_pop = np.real(af.population_from_wf_vector(self.chromosomes[0].prop_psi.propagator_terms.mol.wf.ci))
        return final_pop


    def get_pop_t(self):
        self.chromosomes[0].prop_psi.propagator_terms.mol.wf.set_wf(self.initial_c0, 1)
        self.oc_iterator_parameters.psi_coeff_t = self.chromosomes[0].prop_psi.propagate_n_step(self.oc_iterator_parameters.nstep,
                                                                                                self.chromosomes[0].field.field)
        pop_t = np.real(af.population_from_wf_matrix(self.oc_iterator_parameters.psi_coeff_t))
        return pop_t


    def get_field_t(self):
        return self.oc_iterator_parameters.field_psi_matrix


    def get_field_ampl(self):
        return self.chromosomes[0].field.parameters['fi']

    def get_restart(self):
        out = np.concatenate((self.omegas_matrix, self.chromosomes[0].field.parameters['fi']))
        return out
























    #faccio esattamente la stessa cosa che per DEAP ma è palese
    def init_general_chromosome(self, molecule, starting_field, pcm):
        chromosome = Chromosome()
        chromosome.init_chromosome(self.oc_iterator_parameters.dt, molecule, starting_field, pcm)
        self.chromosomes = self.n_chromosomes * chromosome
        for i in range(self.n_chromosomes):
            rand = []
            for j in range(self.n_amplitudes):
                rand.append(random.uniform(-0.01, 0.01))
            self.chromosomes[i].amplitudes = deepcopy(rand)



    def evaluate_J_general_chromosomes(self):
        for i in range(len(self.chromosomes)):
            self.chromosomes[i].amplitudes_to_field()
            self.chromosomes[i].prop_psi.propagate_n_step(self.oc_iterator_parameters.nstep, self.chromosomes[i].field.field)
            self.n_mutation_succes += self.chromosomes[i].calc_J(self.oc_iterator_parameters.target_state,
                                                                 self.oc_iterator_parameters.alpha_t,
                                                                 self.oc_iterator_parameters.dt)

#    def evolve_general(self):
#        n_children_each = self.n_chromosomes / self.n_selected_chr
#        if n_children_each.is_integer() == False:
#            sys.exit("wrong numer of chromosome to evolve, n_chromosomes/n_selected must be integer")
#        print("evolve")
#        # Select the next generation individuals
#        selected = self.evolutionary_algorithms.select(self.chromosomes, fit_attr='J')
#        selected = self.evolutionary_algorithms.clone(selected)
#        new = []
#        for i in range(int(n_children_each - n_mutate)):
            # se tolgo questo il primo abbinamento me lo fa tra i vicini di J, se lo lascio è random
#            random.shuffle(selected)
#            for first, second in zip(selected[::2], selected[1::2]):
#                if random.random() < mate_probability:
#                    a, b = self.evolutionary_algorithms.mate(deepcopy(first.amplitudes), deepcopy(second.amplitudes))
#                else:
#                    a, b = deepcopy(first.amplitudes), deepcopy(second.amplitudes)
#                new.append(deepcopy(a))
#                new.append(deepcopy(b))
#            random.shuffle(selected)





