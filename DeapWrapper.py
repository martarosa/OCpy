import random

from deap import base
from deap import creator
from deap import tools

import numpy as np

import PropagatorsEulero as prop
import auxiliary_functions as af


from OCIterator import OCIterator
from Field import Field


folder = "/home/mana/programmi/python/optimal_control/OCpy/test/oc/fullpy/"
namefile = "input.dat"


OC_system=ini.SystemManager()
OC_system.init_system(folder, namefile)




class Chromosome():
    def __init__(self):
        self.ampl_genes = []  # n_ampitudes
        self.J = None

        self.field = Field()
        self.prop_psi = prop.PropagatorEulero2Order()


    def init_chromosome(self, dt, molecule, field, pcm):
        self.field = field
        self.prop_psi.set_propagator(dt, molecule, pcm)
        self.ampl_genes = np.concatenate((self.field.parameters['fi'], self.field.parameters['fi_cos']), axis = 1)








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
