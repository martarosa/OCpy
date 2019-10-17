import SystemManager as ini
import argparse
import time as time
from deap import base, creator
import numpy as np


import random

from deap import base
from deap import creator
from deap import tools




folder = "/home/mana/programmi/python/optimal_control/OCpy/test/2order_prop/vac/26-9/"
namefile = "input_nuovo.dat"


OC_system=ini.SystemManager()
OC_system.init_system(folder, namefile)

#creo un oggetto di tipo fitness, che è implementato. da solo non fa niente, è un attributo assegnato a un individuo

creator.create("FitnessMax", base.Fitness, weights=(1.0,))

#creo il mio cromosoma. è una lista di numeri.

creator.create("Individual", list, fitness=creator.FitnessMax)
# questa è una classe che si comporta come una lista e ha un atrtibuto fitness di tipo FitnessMax
# posso anche farlo comportare come un array numpy o un array.array
# è come dire   class Individual()
#                  def __init__(self):
#                      self.fitness = FitnessMax()
#e poi, visto che è uan lista, una cosa del tipo
#                      self.nome = []
#che fa sì che

i = creator.Individual()
print(i)
print("mi ha stampato una lista vuota")
i.append([1,2,3])
print(i)
print("me l'ha riempita")
i.fitness.values=[40]
print(i.fitness.values)
print("mi ha riempito anche la fitness")

size = 10

#per fare le funzioni invece non usiamo creator che va bene per le classi ma toolbox per le funzioni
toolbox = base.Toolbox()
toolbox.register("printTB", print)
toolbox.printTB("mana")
#mi stampa mana perchè manda alla funzione print
# è come def printTB(name):
#               print(name)
toolbox.register("attr_float", random.random)
#è come chiamare random
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, n=size)
# questa è una funzione individual, che funziona come initRepeat e le dà come argomenti i tre successivi.
#initrepeat si aspetta 3 argomenti: un container in cui mettere quello che crea, che vuol dire la cosa che ritorna
#                                   una funzione con cui creare
#                                   quante volte chiamare questa funzione

#quindi se ora dico
a = toolbox.individual()
print(a)
print("a e un Individual riempito")
print(a.fitness.values)
print("senza valori di fitness")

##############
#Permutation
#############

toolbox.register("indices_new_order", range(size), size)
#random.sample prende una lista e un intero k e restituisce una lista di lughezza k con elementi prsi dalla lista data
# quindi qui la lista è range(10) [0,1,2,3,4,5,6,7,8,9] e ne sceglie random 10 di questi SENZA ripetizioni quindi
# indices_new_order mi sta dando un nuovo ordine degli indici. Che mi serve per fare la permutazione!!!
toolbox.indices_new_order()

toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.indices_new_order)
#initIterate prende due argomenti. Ritorna il tipo del primo (Individual()) e lo fa con la funzione secondo argomento
#quindi rispetto alla riga sopra invece di fare una lista e stop la caccia in un tipo Iterator

##########################
#EvolutionStrategies
##########################










