from OC.OCEuleroIterator import Eulero1PropagationIterator, Eulero2PropagationIterator
from OC.OCGeneticIterator import OCGeneticIterator
from OC.OCRabitzIterator import OCRabitzIterator
from medium.VacMedium import VacMedium
from medium.FrozenSolventMedium import FrozenSolventMedium
from medium.DinamicMedium import DinamicMedium
from propagator.ClassicalPropagatorTerms import ClassicalPropagatorTerms
from propagator.PropagatorsEulero import PropagatorEulero2Order, PropagatorEulero1Order
from propagator.PropagatorOCRabitz import PropagatorOCbwd, PropagatorOCfwd
from save.SaveEulero import SaveEulero
from save.SaveOCRabitz import SaveOCRabitz
from save.SaveOCGenetic import SaveOCGenetic



OCAlgorithmDict = {
                    "genetic": OCGeneticIterator,
                    "rabitzi": OCRabitzIterator,
                    "rabitzii": OCRabitzIterator,
                    "eulero_2order_prop": Eulero2PropagationIterator,
                    "eulero_1order_prop": Eulero1PropagationIterator
}


SaveDict = {
    "genetic": SaveOCGenetic,
    "rabitzi": SaveOCRabitz,
    "rabitzii":SaveOCRabitz ,
    "eulero_2order_prop": SaveEulero,
    "eulero_1order_prop": SaveEulero
}


MediumDict = {
    "vac": VacMedium,
    "sol": FrozenSolventMedium,
    "nanop": DinamicMedium
}

PropagatorDict = {
    "eulero_2order_prop": PropagatorEulero2Order,
    "eulero_1order_prop": PropagatorEulero1Order,
    "rabitzi": PropagatorOCfwd,
    "rabitzii": PropagatorOCfwd
}

PropagatorTermsDict = {
    "classical": ClassicalPropagatorTerms,

}

