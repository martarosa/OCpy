import numpy as np

from read_and_set.read import auxiliary_functions as af
from OC.ABCOCIterator import ABCOCIterator, OCIteratorParameters
from SystemObj import DiscreteTimePar
from Chromosome import Chromosome
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


class GeneticParameters():
    def __init__(self):
        self.n_chromosomes = None
        self.n_selected_chr = None

        self.amplitude_lim = None
        self.n_amplitudes = None
        self.omegas_matrix = None
        self.mate = None
        self.mate_probability = None

        self.mutate = None
        self.mutate_probability = None
        self.n_mutate = None
        self.mutate_mu = None
        self.mutate_starting_sigma = None
        self.eta_thr = None
        self.q = None
        self.select = None


class OCGeneticIterator(ABCOCIterator):
    def __init__(self):
        super().__init__()
        self.par = OCIteratorParameters()
        self.simulation_par = DiscreteTimePar()
        self.genetic_par = GeneticParameters()

        self.convergence_t = None
        self.J = None
        self.field_psi_matrix = Func_tMatrix()
        self.psi_coeff_t_matrix = None
        self.dict_out = {}

        #self.prop_psi = prop.PropagatorEulero2Order()

        #genetic
        self.initial_c0 = None
        self.chromosomes = [] #can be a DEAP class or personal chromosome class




    def iterate(self, current_iteration):
        print(current_iteration)
        self.evolve()
        self.check_convergence()

    def check_convergence(self):
        J_prev_tmp = np.copy(self.J)
        self.calc_J()
        self.convergence_t = self.J - J_prev_tmp
        self.convergence_t = self.convergence_t[0]

    def calc_J(self):
        self.chromosomes.sort(key=lambda x: x.J, reverse=True)
        self.J = self.chromosomes[0].J.values
        self.field_psi_matrix = self.chromosomes[0].field.field


    def init_output_dictionary(self):
        self.dict_out['log_file'] = self.get_log_file_out
        self.dict_out['final_pop'] = self.get_final_pop
        self.dict_out['pop_t'] = self.get_pop_t
        self.dict_out['field_t'] = self.get_field_t
        self.dict_out['field_ampl'] = self.get_field_ampl



    def init(self, molecule, starting_field, pcm, alpha_t, oc_input, iterator_config_input):
        self.simulation_par.nstep = oc_input.nstep
        self.simulation_par.dt = oc_input.dt
        self.par.target_state = oc_input.target_state
        self.par.alpha_t = alpha_t

        self.convergence_t = 99999
        self.J = 99999

        self.field_psi_matrix = deepcopy(starting_field.field)

        self.psi_coeff_t_matrix = np.zeros([self.simulation_par.nstep + 1, molecule.wf.n_ci], dtype=complex)
        self.init_output_dictionary()

#        self.psi_coeff_t = np.zeros([self.simulation_par.nstep + 1, self.prop_psi.propagator_terms.mol.wf.n_ci + 1], dtype=complex)
#        self.psi_coeff_t[:, 0] = np.linspace(0, self.simulation_par.dt * self.simulation_par.nstep , self.simulation_par.nstep + 1)
#        self.psi_coeff_t[0, 1:] =  self.initial_c0

        self.init_genetic(molecule, starting_field, pcm, iterator_config_input)


    def init_genetic(self, molecule, starting_field, pcm, genetic_input):
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

        self.init_chromosomes(molecule, starting_field, pcm)
    #qui inizializzo l'algoritmo genetico che sto usando, che per adesso è solo uno
        self.init_evolutionary_algorithm(genetic_input)



    def init_chromosomes(self, molecule, starting_field, pcm):
            self.init_DEAP_chromosomes(molecule, starting_field, pcm)

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
        prop_psi.set_propagator(molecule, pcm)
        self.initial_c0 = prop_psi.propagator_terms.mol.wf.ci
    #... e dentro a ogni cromosoma della mia lista inizializzo il campo, il propagatore e metto J = 0.5
    #Sono tutti uguali ma le ampiezze del mio vettore sono diverse perchè le ho fatte random prima
        for i in range(len(self.chromosomes)):
            self.chromosomes[i].field = deepcopy(starting_field)
            self.chromosomes[i].prop_psi = deepcopy(prop_psi)
            self.chromosomes[i].J.values = [0.5]
        self.par.J = [99999]


    def evaluate_J_DEAP_chromosome(self, chro):
    # prende l'array di Chromosme e lo mette nel campo come ampiezze
        chro.field.par.fi = np.asarray(chro).reshape((-1, 3))
    #genera il campo con forma sum
        chro.field.chose_field('sum_pip', simulation_par = self.simulation_par)
        chro.field.chose_field()
        #chro.field.chose_field('sum_pip')
    #setto il propagatore e lo propago con il campo appena generato
        chro.prop_psi.propagator_terms.mol.wf.set_wf(self.initial_c0, 1)
        chro.prop_psi.propagate_n_step(self.simulation_par.nstep, chro.field.get_full_field())
    #calcolo popolazione e integrale del campo, per debug
        pop= np.real(af.projector_mean_value(chro.prop_psi.propagator_terms.mol.wf.ci,
                                             self.par.target_state))
        field = np.real(self.alpha_field_J_integral(chro.field,
                                                  self.par.alpha_t))
    #calcolo J
        J = np.real(af.projector_mean_value(chro.prop_psi.propagator_terms.mol.wf.ci,
                                            self.par.target_state)
                    - self.alpha_field_J_integral(chro.field,
                                                self.par.alpha_t))

        success =  J - chro.J.values[0]
        if success > 0:
            self.n_mutation_succes +=1
        return [J]


    def evolve(self):
        self.evolve_subsequent_DEAP()


    def evolve_subsequent_DEAP(self):
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
        self.check_bounds(new)
        J = self.genetic_algorithms.map(self.genetic_algorithms.evaluate, new)
        for ind, fit in zip(new, J):
            ind.J.values = fit
        self.chromosomes[:] = new


    #quello del paper, con parte crossover e parte mutati
    def evolve_mixed_DEAP(self):
        n_mate = self.genetic_par.n_chromosomes-(self.genetic_par.n_selected_chr * self.genetic_par.n_mutate)
        print("evolve")
        # Select the next generation individuals
        new_crossover = []
    #l'algoritmo è già definito. Quindi seleziono da cromosomi i migliori secondo J
        selected = self.genetic_algorithms.select(self.chromosomes, fit_attr='J')
        selected = self.genetic_algorithms.clone(selected)
    # Apply crossover on the offspring
    # li prendo tutti e ne genero quanti me ne servono mescolandoli
        for i in range(n_mate):
    #se tolgo questo il primo abbinamento me lo fa tra i vicini di J, se lo lascio è random
            random.shuffle(selected)
            for first, second in zip(selected[::2], selected[1::2]):
                if random.random() < self.genetic_par.mate_probability:
                    a, b = self.genetic_algorithms.mate(deepcopy(first), deepcopy(second))
                else:
                    a, b = deepcopy(first), deepcopy(second)
                new_crossover.append(deepcopy(a))
                new_crossover.append(deepcopy(b))
    # Apply mutation on the offspring
    # li prendo tutti in ordine e li muto le volte che mi serve
        new_mutation = []
        print("sigma:")
        print(self.genetic_algorithms.mutate.keywords['sigma'])
        for i in range(self.genetic_par.n_selected_chr):
            for j in range(self.genetic_par.n_mutate):
                new_mutation.append(deepcopy(selected[i]))
                self.genetic_algorithms.mutate(new_mutation[-1])
        #self.check_bounds(new_mutation)

        #calcolo la J del crossover e la J della mutazione
        J_crossover = list(self.genetic_algorithms.map(self.genetic_algorithms.evaluate, new_crossover))
        self.n_mutation_succes = 0
        J_mutation = list(self.genetic_algorithms.map(self.genetic_algorithms.evaluate, new_mutation))
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
        eta = self.n_mutation_succes/(self.genetic_par.n_mutate * self.genetic_par.n_selected_chr)
        print(eta)
        #print(eta)
        #print(self.mutate_sigma)
        if eta <= self.genetic_par.eta_thr:
            self.mutate_sigma = self.mutate_sigma*q
        else:
            self.mutate_sigma = self.mutate_sigma / q
        self.genetic_algorithms.mutate.keywords['sigma'] = self.mutate_sigma
        self.n_mutation_succes = 0

    def check_bounds(self, matrix):
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                if matrix[i][j] > self.genetic_par.amplitude_lim:
                    matrix[i][j] = self.genetic_par.amplitude_lim
                elif matrix[i][j] < -self.genetic_par.amplitude_lim:
                    matrix[i][j] = -self.genetic_par.amplitude_lim
                #elif (matrix[i][j] < 0.0005 and matrix[i][j] > -0.0005):
                #    matrix[i][j] = 0.0

    def init_evolutionary_algorithm(self, iterator_parameters):
        self.genetic_algorithms = base.Toolbox()
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
        self.genetic_algorithms.register('evaluate',
                                         self.evaluate_J_DEAP_chromosome)




    def get_log_file_out(self):
        integral_field = np.real(self.field_J_integral(self.field_psi_matrix, self.simulation_par.dt))
        norm_proj = np.real(af.projector_mean_value(self.chromosomes[0].prop_psi.propagator_terms.mol.wf.ci,
                                                    self.par.target_state) /
                            (np.dot(self.chromosomes[0].prop_psi.propagator_terms.mol.wf.ci,
                                    np.conj(self.chromosomes[0].prop_psi.propagator_terms.mol.wf.ci))))
        return np.array([self.convergence_t, self.J[0], norm_proj, integral_field])


    def get_final_pop(self):
        final_pop = np.real(af.population_from_wf_vector(self.chromosomes[0].prop_psi.propagator_terms.mol.wf.ci))
        return final_pop


    def get_pop_t(self):
        self.chromosomes[0].prop_psi.propagator_terms.mol.wf.set_wf(self.initial_c0, 1)
        self.par.psi_coeff_t = self.chromosomes[0].prop_psi.propagate_n_step(self.simulation_par.nstep,
                                                                             self.chromosomes[0].field.field)
        pop_t = np.real(af.population_from_wf_matrix(self.par.psi_coeff_t))
        return pop_t


    def get_field_t(self):
        return self.field_psi_matrix


    def get_field_ampl(self):
        return self.chromosomes[0].field.parameters['fi']

    def get_restart(self):
        out = np.concatenate((self.genetic_par.omegas_matrix, self.chromosomes[0].field.parameters['fi']))
        return out



    #faccio esattamente la stessa cosa che per DEAP ma è palese
    def init_general_chromosome(self, molecule, starting_field, pcm):
        chromosome = Chromosome()
        chromosome.init_chromosome(self.par.dt, molecule, starting_field, pcm)
        self.chromosomes = self.genetic_par.n_chromosomes * chromosome
        for i in range(self.genetic_par.n_chromosomes):
            rand = []
            for j in range(self.genetic_par.n_amplitudes):
                rand.append(random.uniform(-0.01, 0.01))
            self.chromosomes[i].amplitudes = deepcopy(rand)



    def evaluate_J_general_chromosomes(self):
        for i in range(len(self.chromosomes)):
            self.chromosomes[i].amplitudes_to_field()
            self.chromosomes[i].prop_psi.propagate_n_step(self.par.nstep, self.chromosomes[i].field.field)
            self.n_mutation_succes += self.chromosomes[i].calc_J(self.par.target_state,
                                                                 self.par.alpha_t,
                                                                 self.par.dt)

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





