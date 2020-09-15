from propagator import PropagatorsEulero as prop
from propagator import PropagatorQuantum as qprop




PropagatorDictionary ={ "eulero_1order": prop.PropagatorEulero1Order(),
                        "eulero_2order": prop.PropagatorEulero2Order(),
                        "quantum":       qprop.PropagatorQuantum()}