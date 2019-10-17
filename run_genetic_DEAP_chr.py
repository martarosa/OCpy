import numpy as np
from copy import deepcopy

import random

from deap import base
from deap import creator
from deap import tools
from deap import algorithms



###################DEAP INITIALIZATION#################
#mate
mate_probability = 1 #given 2 list which is the probability they mate
#cxUniform
cx_uniform_probability = 0.5
#mutate
mutate_probability = 0.1 # applied to each element of array
#mutGaussian
mu = 0.1
sigma = 0.1
#######################################################



target_state = np.array([0,1])
n_chromosomes = 48
n_selected_chr = 8
n_to_generate_from_each = int(n_chromosomes/n_selected_chr) + 1

n_amplitudes = 4
ampl_min = -2
ampl_max = 2


def create_deap_chromosomes():
    creator.create("J", base.Fitness, weights=(-1.0,))
    creator.create("Chromosome", list, fitness=creator.J, prop_psi=[])
    toolbox = base.Toolbox()
    toolbox.register("create_random_ampl", random.uniform, ampl_min, ampl_max)
    toolbox.register("single_chromosome", tools.initRepeat, creator.Chromosome, toolbox.create_random_ampl,
                          n=n_amplitudes)
    toolbox.register('chromosomes_population', tools.initRepeat, list, toolbox.single_chromosome)
    chromosomes = toolbox.chromosomes_population(n=n_chromosomes)  # create all chromosomes
    return chromosomes


def evaluate_DEAP_chromosomes(chro):
    for i in range(len(chro)):
        matr = np.asarray(chro[i]).reshape(2, 2)
        chro[i].fitness.values = diff(matr)
    chro.sort(key=lambda x: x.fitness.values, reverse=True)


def evaluate_single(chro):
    matr = np.asarray(chro).reshape(2, 2)
    return diff(matr)


def diff(matr):
    zz = np.sqrt((0-matr[0][0])*(0-matr[0][0]))
    zu = np.sqrt((1 - matr[0][1]) * (1 - matr[0][1]))
    uz = np.sqrt((1 - matr[1][0]) * (1 - matr[1][0]))
    uu = np.sqrt((0 - matr[1][1]) * (0 - matr[1][1]))
    return [zz+zu+uz+uu]

def evolve_DEAP(chro):
    # Select the next generation individuals
    selected = toolbox.select(chro)
    new = []
    # Clone the selected individuals
    selected = toolbox.clone(selected)
    # Apply crossover on the offspring
    for i in range(n_to_generate_from_each):
        for first, second in zip(selected[::2], selected[1::2]):
            if random.random() < mate_probability:
                a, b = toolbox.mate(deepcopy(first), deepcopy(second))
                new.append(deepcopy(a))
                new.append(deepcopy(b))
        random.shuffle(selected)
    new = new[:n_chromosomes]

    # Apply mutation on the offspring
    for mutant in new:
        if random.random() < mutate_probability:
            toolbox.mutate(mutant)
    check_bounds(new)
    return new


def total_EA(chro):
    new = []
    selected = toolbox.select(chro)
    selected = toolbox.clone(selected)
    # Apply crossover on the offspring
    for i in range(n_to_generate_from_each):
        for first, second in zip(selected[::2], selected[1::2]):
            if random.random() < mate_probability:
                a, b = toolbox.mate(deepcopy(first), deepcopy(second))
                new.append(deepcopy(a))
                new.append(deepcopy(b))
        random.shuffle(selected)
    new = new[:n_chromosomes]

    # Apply mutation on the offspring
    for mutant in new:
        if random.random() < mutate_probability:
            toolbox.mutate(mutant)
    check_bounds(new)
    fitness = toolbox.map(toolbox.evaluate, new)
    for ind, fit in zip(new, fitness):
        ind.fitness.values = fit
    return new


def check_bounds(chro):
    for i in range(len(chro)):
        for j in range(n_amplitudes):
            if chro[i][j] > ampl_max:
                chro[i][j] = ampl_max
            elif chro[i][j] < ampl_min:
                chro[i][j] = ampl_min

def init_DEAP_evolutionary_algorithm():
    toolbox.register("mate", tools.cxUniform, indpb=cx_uniform_probability)
    toolbox.register("mutate", tools.mutGaussian, mu = mu, sigma = sigma, indpb = mutate_probability)
    toolbox.register("select", tools.selBest, k = n_selected_chr)
    toolbox.register('evaluate', evaluate_single)

toolbox = base.Toolbox()
init_DEAP_evolutionary_algorithm()

chromosomes = create_deap_chromosomes()
evaluate_DEAP_chromosomes(chromosomes)




for i in range(1000):
    chromosomes[:] = total_EA(chromosomes)



#for i in range(1000):
#    tmp = evolve_DEAP(chromosomes)
#    for i in range(len(chromosomes)):
#        chromosomes[i] = deepcopy(tmp[i])
#    evaluate_DEAP_chromosomes(chromosomes)



















