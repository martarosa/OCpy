class MoleculeParameters:
    def __init__(self):
        self.muT = None  # transition dipoles from external code
        self.en_ci = None  # excitation energies from external code
        self.Vijn = None  # potential from external code if PCM

        self.Vijn_fortran_flip =  None