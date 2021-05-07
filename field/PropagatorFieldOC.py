import numpy as np
from read_and_set.read import auxiliary_functions as af



#used for Rabitz optimization algorithm
#Field (x y z) components at time t: 3D vector
#return 1 x 3 vector : field at time t+1

class PropagatorFieldOC():
    def __init__(self):
        self.field_dt_vector = np.array(3)

    def propagate_field_OC_Rabitz(self, ket, bra, muT, alpha):
        f = af.double_summation(ket, np.conj(bra), muT)
        self.field_dt_vector = -1 * np.imag(f) / alpha
        return self.field_dt_vector
