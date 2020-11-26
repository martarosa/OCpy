from propagator.PropagatorOCRabitz import PropagatorOCfwd
from propagator.PropagatorsEulero import PropagatorEulero2Order, PropagatorEulero1Order, PropagatorEulero2OrderGeneralPerturbation
from propagator.QuantumPropagatorTerms import QuantumPropagatorTerms
from propagator.ClassicalPropagatorTerms import ClassicalPropagatorTerms
from propagator.QuantumPropagator import PropagatorQuantum
from read_and_set.read.conf.ReadNoConf import ReadNoConf
from read_and_set.set.SetNoConf import SetNoConf
from read_and_set.read.conf.ReadIBMConf import ReadIBMConf
from read_and_set.set.SetIBMInput import SetIBMInput


PropagatorDict = {
    "eulero_2order": PropagatorEulero2Order,
    "eulero_1order": PropagatorEulero1Order,
    "eulero_2order_psi4": PropagatorEulero2OrderGeneralPerturbation,
    "rabitz": PropagatorOCfwd,
    "quantum_trotter_suzuki" : PropagatorQuantum
    
}
PropagatorTermsDict = {
    "eulero_2order": ClassicalPropagatorTerms,
    "eulero_1order": ClassicalPropagatorTerms,
    "rabitz": ClassicalPropagatorTerms,
    "eulero_2order_psi4": ClassicalPropagatorTerms,
    "quantum_trotter_suzuki" : QuantumPropagatorTerms
}

PropagatorConfig = {
        "eulero_2order" : ReadNoConf,
        "eulero_1order" : ReadNoConf,
        "rabitz" : ReadNoConf,
        "eulero_2order_psi4": ReadNoConf,
        "quantum_trotter_suzuki" : ReadIBMConf
}

PropagatorSet = {
        "eulero_2order" : SetNoConf,
        "eulero_1order" : SetNoConf,
        "rabitz" : SetNoConf,
        "eulero_2order_psi4": SetNoConf,
        "quantum_trotter_suzuki" : SetIBMInput
}