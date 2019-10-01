import numpy as np


class MolecularHamiltonian:
    def __init__(self):
        self.muT = None  # transition dipoles from external code
        self.en_ci = None  # excitation energies from external code
        self.toSubtract = None

    def set_muT(self, muT):
        self.muT = muT

    def mod_muT_local_field(self, mu_local_field):
        self.muT += np.real(mu_local_field)


    def set_en_ci(self, en_ci):
        self.en_ci = en_ci

    def get_muT(self):
        return self.muT

    def get_en_ci(self):
        return self.en_ci


