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




class Chromosome():
    def __init__(self):
        self.amplitudes = []  # n_ampitudes
        self.J = None
        self.prop_psi = []

def create_chromosomes():
    chromosomes = []
    starting_chr = Chromosome()
    starting_chr.prop_psi = [1, 0]

    for i in range(n_chromosomes):
        chromosomes.append(deepcopy(starting_chr))
        rand = []
        for j in range(n_amplitudes):
            rand.append(random.uniform(ampl_min, ampl_max))
        chromosomes[-1].amplitudes = deepcopy(rand)
    return chromosomes

def evaluate_chromosomes(chro):
    for i in range(len(chro)):
        matr = np.asarray(chro[i].amplitudes).reshape(2, 2)
        chro[i].J = diff(matr)
    chro.sort(key=lambda x: x.J, reverse=True)

def diff(matr):
    zz = np.sqrt((0-matr[0][0])*(0-matr[0][0]))
    zu = np.sqrt((1 - matr[0][1]) * (1 - matr[0][1]))
    uz = np.sqrt((1 - matr[1][0]) * (1 - matr[1][0]))
    uu = np.sqrt((0 - matr[1][1]) * (0 - matr[1][1]))
    return zz+zu+uz+uu

def evolve(chro):
    # Select the next generation individuals
    selected = toolbox.select(chro,  fit_attr='J')
    new = []
    # Clone the selected individuals
    selected = toolbox.clone(selected)
    # Apply crossover on the offspring
    for i in range(n_to_generate_from_each):
        for first, second in zip(selected[::2], selected[1::2]):
            if random.random() < mate_probability:
                a, b = toolbox.mate(deepcopy(first.amplitudes), deepcopy(second.amplitudes))
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
    toolbox.register("select", tools.selWorst, k = n_selected_chr)
    toolbox.register('evaluate', evaluate_chromosomes)

toolbox = base.Toolbox()
init_DEAP_evolutionary_algorithm()

chromosomes = create_chromosomes()
evaluate_chromosomes(chromosomes)

for i in range(4000):
    tmp = evolve(chromosomes)
    for i in range(len(chromosomes)):
        chromosomes[i].amplitudes = deepcopy(tmp[i])
    evaluate_chromosomes(chromosomes)










