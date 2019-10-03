from deap import base
from deap import creator
from deap import tools

import numpy as np

from propagator import PropagatorsEulero as prop

from field.Field import Field
from Chromosome import Chromosome

import random




def deap_creator():

     creator.create("J", base.Fitness, weights=(1.0,))
     creator.create("Chromosome", list, fitness = creator.J)
     toolbox = base.Toolbox()
     toolbox.register("create_random_ampl", random.uniform, -0.001, 0.001)
     toolbox.register("single_chromosome", tools.initIterate, creator.Chromosome, toolbox.create_random_ampl, n=1)





#class EvolutionaryStrategy():
#    def __init__(self):
#        self.evolutionary_algorithms = EvolutionaryOperators()
#        self.evolutionary_strategy = []




#class EvolutionaryOperators():
#    def __init__(self):

#        self.dict_terms = {}


#    def init_terms_dictionary(self):
#        self.dict_terms["1pt_crossover"] = self.eulero1_coeff_term
#        self.dict_terms["2pt_crossover"] = self.eulero2_coeff_term
#        self.dict_terms["uniform_crossover"] = self.eulero_energy_term
#        self.dict_terms["gaussian_mutation"] = self.eulero_field_term
