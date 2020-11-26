import numpy as np

class MoleculeInput:
    def __init__(self):
        self.wf_ci = None
        self.muT = None
        self.en_ci = None
        self.Vijn = np.array([None])
        self.nmo = None
        self.ao_coefficients = None
        self.two_electron_int = None
        self.ao_kinetic_integral = None
        self.ao_potential_integral = None