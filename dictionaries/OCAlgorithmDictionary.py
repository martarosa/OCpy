
from read_and_set.read.ReadNamelistGenetic import ReadNamelistGenetic

from read_and_set.set.SetGeneticOCInput import SetGeneticOCInput
from read_and_set.set.SetNoConfOCInput import SetNoConfOCInput

from deap import tools



OCAlgorithmDictionary ={"genetic": ReadNamelistGenetic,
                        "rabitzi": None,
                        "rabitzii":None}


SetOCInputDictionary = {"genetic":   SetGeneticOCInput,
                        "rabitzi":   SetNoConfOCInput,
                        "rabitzii":  SetNoConfOCInput}


EvolutionaryAlgorithmDictionary = {'cxUniform':   tools.cxUniform,
                                   'mutGaussian': tools.mutGaussian,
                                   'selBest':     tools.selBest}