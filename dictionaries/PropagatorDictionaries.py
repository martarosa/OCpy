from propagator.PropagatorOCRabitz import PropagatorOCfwd
from propagator.PropagatorsEulero import PropagatorEulero2Order, PropagatorEulero1Order

PropagatorDict = {
    "eulero_2order": PropagatorEulero2Order,
    "eulero_1order": PropagatorEulero1Order,
    "rabitz": PropagatorOCfwd,
}


