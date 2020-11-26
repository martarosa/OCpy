from propagator.ClassicalPropagatorTerms import ClassicalPropagatorTerms
from propagator.QuantumPropagatorTerms import QuantumPropagatorTerms

PropagatorTermsDict = {
    "eulero_2order": ClassicalPropagatorTerms,
    "eulero_1order": ClassicalPropagatorTerms,
    "rabitz": ClassicalPropagatorTerms,
    "eulero_2order_psi4": ClassicalPropagatorTerms,
    "quantum_trotter_suzuki" : QuantumPropagatorTerms
}