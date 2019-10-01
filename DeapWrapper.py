from deap import base
from deap import creator
from deap import tools

import numpy as np

from propagator import PropagatorsEulero as prop

from field.Field import Field
from Chromosome import Chromosome






class DeapWrapper():


    def init(self, size):
        creator.create("J", base.Fitness, weights=(1.0,))
        creator.create("Chromosome", list, fitness = creator.J, field = Field(), prop_psi = prop.PropagatorEulero2Order())

        toolbox = base.Toolbox()
        #toolbox.register("attr_float", random.random)


        toolbox.register("individual", tools.initRepeat, creator.Individual,
                         toolbox.attr_float, n=IND_SIZE)



class EvolutionaryStrategy():
    def __init__(self):
        self.evolutionary_algorithms = EvolutionaryOperators()
        self.evolutionary_strategy = []




class EvolutionaryOperators():
    def __init__(self):

        self.dict_terms = {}


    def init_terms_dictionary(self):
        self.dict_terms["1pt_crossover"] = self.eulero1_coeff_term
        self.dict_terms["2pt_crossover"] = self.eulero2_coeff_term
        self.dict_terms["uniform_crossover"] = self.eulero_energy_term
        self.dict_terms["gaussian_mutation"] = self.eulero_field_term
