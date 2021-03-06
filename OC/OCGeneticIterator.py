import time
import itertools as itt
import numpy as np
import pandas as pd
import multiprocessing as mp
import concurrent.futures
import random

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

    def iterate(self, current_iteration):
        print(current_iteration)
        self.evolve()
        self.check_convergence()

    def init_output_dictionary(self):
        self.dict_out['log_file'] = self.get_log_file_out
        self.dict_out['final_pop'] = self.get_final_pop
        self.dict_out['pop_t'] = self.get_pop_t
        self.dict_out['field_t'] = self.get_field_t
        self.dict_out['field_ampl'] = self.get_field_ampl

    def init(self, molecule, starting_field, starting_medium, alpha_t, oc_input, iterator_config_input):
        self.discrete_t_par.nstep = oc_input.nstep
        self.discrete_t_par.dt = oc_input.dt
        self.par.target_state = oc_input.target_state
        self.par.alpha_t = alpha_t
        self.par.convergence_t = 99999
        self.par.J = 99999
        self.medium = starting_medium
        self.molecule = molecule
        self.start_field = starting_field
        self.field_psi_matrix = deepcopy(starting_field.field)
        self.init_output_dictionary()
        self.init_genetic(molecule, starting_field, starting_medium, iterator_config_input)

    def init_genetic(self, molecule, starting_field, starting_medium, genetic_input):
        self.genetic_par.genetic_algorithm = genetic_input.genetic_algorithm
        self.genetic_par.n_chromosomes = genetic_input.n_chromosomes
        self.genetic_par.n_selected_chr = genetic_input.n_selected_chr
        self.genetic_par.amplitude_lim = genetic_input.amplitude_lim

        self.genetic_par.mate = genetic_input.mate
        self.genetic_par.mate_probability = genetic_input.mate_probability

        self.genetic_par.mutate = genetic_input.mutate
        self.genetic_par.mutate_probability = genetic_input.mutate_probability
        self.genetic_par.mutate_mu = genetic_input.mutate_mu
        self.genetic_par.mutate_starting_sigma = genetic_input.mutate_starting_sigma
        self.genetic_par.select = genetic_input.select
        self.genetic_par.n_amplitudes = starting_field.par.fi.size
        self.genetic_par.omegas_matrix = starting_field.par.omega
        self.init_chromosomes(molecule, starting_field, starting_medium)
        self.init_evolutionary_algorithm(genetic_input)

    def init_chromosomes(self, molecule, starting_field, starting_medium):
        print("init_chromo")
        toolbox= base.Toolbox()
    # fitness is an abstract class.
    # J is a child of fitness with positive weight, which means that whants to be maximized
        creator.create("J", base.Fitness, weights=(1.0,))
    # here we create an object called Chromosome, which has as attributes J, a Field() object and a Propagator() object.
    # Chromosome() objact has also an empty list associated. When you mate, mutate, select the cromosome with deap library is the list that it is worked on.
        creator.create("Chromosome", list, J=creator.J, field = Field(), prop_psi = prop.PropagatorEulero2Order())
    #here I create a single instance of Chormosome() type.
    # All the values in the list are set to 0
        toolbox.register("single_chromosome",
                         tools.initRepeat, creator.Chromosome,
                         self.init_chromosomes_zero,
                         n=self.genetic_par.n_amplitudes)
    #the we create a list of chromosomes all zero valued
        toolbox.register('chromosomes_population',
                         tools.initRepeat,
                         list,
                         toolbox.single_chromosome)
        self.chromosomes = toolbox.chromosomes_population(n=self.genetic_par.n_chromosomes) #create all chromosomes
        #----------------
        #init the propagation part
        #----------------
        self.initial_c0 = molecule.wf.ci
        # all the Field() objects in the Chromosomes are initialized copying starting_field = Field() object,
        # which provides the omegas. amplitudes are all equal to the default value 0.01 and are then copied in the list
        # and mutated
        # before starting the first optimization iteration
        self.init_chromosomes_values()
        for i in range(len(self.chromosomes)):
            self.chromosomes[i].field = deepcopy(starting_field)
            #ggiorno le ampiezze del campo con i calori generati random del cromosoma
            if (starting_field.par.field_type == "restart_genetic"):
                self.field_amplitudes_to_chromosome(self.chromosomes[i], starting_field.par.fi.reshape((1, -1)))
            self.chromosome_coefficients_to_field(self.chromosomes[i])
            medium = deepcopy(starting_medium)
            medium.reset_medium(molecule, self.chromosomes[i].field.field)
            prop_psi = prop.PropagatorEulero2Order()
            prop_psi.set_propagator(molecule, medium)
            self.chromosomes[i].prop_psi = deepcopy(prop_psi)
            self.chromosomes[i].J.values = [0.5]
        self.par.J = [99999]



    def evolve(self):
        print("evolve")
        n_children_each = self.genetic_par.n_chromosomes / self.genetic_par.n_selected_chr
        if n_children_each.is_integer() == False:
            n_children_each = int(n_children_each) + 1
        else:
            n_children_each = int(n_children_each)
        # Select the next generation individuals
        new = []
        selected = self.genetic_algorithms.select(self.chromosomes, fit_attr='J')
        selected = self.genetic_algorithms.clone(selected)
        # Apply crossover on the offspring
        for _ in itt.repeat(None, n_children_each):
            selected_odd = deepcopy(selected[::2])
            selected_even = deepcopy(selected[1::2])
            if self.genetic_par.n_selected_chr % 2 != 0:
                selected_even.append(random.choice(selected))
            for first, second in zip(selected_odd, selected_even):
                if random.random() < self.genetic_par.mate_probability:
                    a, b = self.genetic_algorithms.mate(first, second)
                else:
                    a, b = first, second
                new.append(deepcopy(a))
                new.append(deepcopy(b))
            np.random.shuffle(selected)
        #for _ in itt.repeat(None,int(self.genetic_par.n_chromosomes/2)+1):
        #    np.random.shuffle(selected)
        ##    if random.random() < self.genetic_par.mate_probability:
        #        a, b = self.genetic_algorithms.mate(deepcopy(selected[0]), selected[1])
        #    else:
        #        a, b = deepcopy(selected[0]), deepcopy(selected[1])
        #    new.append(deepcopy(a))
        #    new.append(deepcopy(b))
        new = new[:self.genetic_par.n_chromosomes]
        # Apply mutation on the offspring
        for mutant in new:
            if random.random() < self.genetic_par.mutate_probability:
                self.genetic_algorithms.mutate(mutant)
                for i in range(len(mutant)):
                    if random.random() < self.genetic_par.mutate_probability:
                        mutant[i]=0.0
        self.check_bounds_matrix(new)

        with concurrent.futures.ProcessPoolExecutor() as executor:
            new = list(executor.map(self.genetic_algorithms.evaluate, new))
        #new = list(self.genetic_algorithms.map(self.genetic_algorithms.evaluate, new))
        self.chromosomes[:] = new


    def reset(self, chro):
        chro.field.par.fi = np.asarray(chro).reshape((-1, 3))
        chro.field.chose_field('sum', dt = self.discrete_t_par.dt)
        chro.prop_psi.mol.wf.set_wf(self.initial_c0, 1)
        chro.prop_psi.medium.reset_medium(chro.prop_psi.mol, chro.field.field)



    def evaluate_J_chromosome(self, chro):
        self.reset(chro)
        chro.prop_psi.propagate_n_step(self.discrete_t_par, chro.field.field)
        J = np.real(af.projector_mean_value(chro.prop_psi.mol.wf.ci,
                                            self.par.target_state)
                    - self.alpha_field_J_integral_chromosome(chro.field.field))
        chro.J.values = [J]
        return chro



    def field_amplitudes_to_chromosome(self, chromosome, reshaped_starting_field_fi):
        n = reshaped_starting_field_fi.size
        for j in range(n):
            chromosome[j] = reshaped_starting_field_fi[0][j] + random.uniform(-self.genetic_par.mutate_starting_sigma/10, -self.genetic_par.mutate_starting_sigma/10)

    def chromosome_coefficients_to_field(self, chro):
        chro.field.par.fi = np.asarray(chro).reshape((-1, 3))
        chro.field.chose_field('sum', dt = self.discrete_t_par.dt)

    def init_chromosomes_values(self):
        random.seed(10)
        for i in range(self.genetic_par.n_chromosomes):
            for j in range(self.genetic_par.n_amplitudes):
                self.chromosomes[i][j] = round(random.uniform(-self.genetic_par.amplitude_lim, self.genetic_par.amplitude_lim),4)

    def init_chromosomes_zero(self):
        return 0.0

    def check_bounds_matrix(self, matrix):
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                if matrix[i][j] > self.genetic_par.amplitude_lim:
                    matrix[i][j] = self.genetic_par.amplitude_lim
                elif matrix[i][j] < -self.genetic_par.amplitude_lim:
                    matrix[i][j] = -self.genetic_par.amplitude_lim

    def check_convergence(self):
        J_prev_tmp = np.copy(self.par.J)
        self.calc_J()
        self.par.convergence_t = self.par.J - J_prev_tmp
        self.par.convergence_t = self.par.convergence_t[0]

    def calc_J(self):
        self.chromosomes.sort(key=lambda x: x.J, reverse=True)
        self.par.J = self.chromosomes[0].J.values
        self.field_psi_matrix = self.chromosomes[0].field.field

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
        self.genetic_algorithms.register('evaluate',
                                         self.evaluate_J_chromosome)

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
        self.reset(self.chromosomes[0])
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

