from propagator.ClassicalPropagatorTerms import ClassicalPropagatorTerms
from propagator.QuantumPropagatorTerms import QuantumPropagatorTerms

PropagatorTermsDict = {
    "eulero_2order": ClassicalPropagatorTerms,
    "eulero_1order": ClassicalPropagatorTerms,
    "rabitz": ClassicalPropagatorTerms,
    "quantum_trotter_suzuki" : QuantumPropagatorTerms
}