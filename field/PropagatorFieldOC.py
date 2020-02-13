import numpy as np
from read_and_set.read import auxiliary_functions as af


#-------------------------------------------
#Field (x y z) components at time t: 3D vector
#-------------------------------------------
#return 1 x 3 vector : field at time t+1

class PropagatorFieldOC():
    def __init__(self):
        self.field_dt_vector = np.array(3)

    def propagate_field_OC_Rabitz(self, ket, bra, muT, alpha):
        f = af.double_summation(ket, np.conj(bra), muT)
        self.field_dt_vector = -1 * np.imag(f) / alpha
        return self.field_dt_vector

    def propagate_field_OC_projector(self, coeff_ket, coeff_bra, muT, alpha):
        f_old = np.dot(np.conj(coeff_bra), np.sum(muT*coeff_ket[np.newaxis, :, np.newaxis], axis=1))
        f = np.dot(np.conj(coeff_ket), coeff_bra) * f_old
        self.field_dt_vector = -1 * np.imag(f) / alpha
        return self.field_dt_vector