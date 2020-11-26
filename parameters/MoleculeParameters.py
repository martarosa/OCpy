import numpy as np

class MoleculeParameters:
    def __init__(self):
        self.muT = None  # transition dipoles from external code
        self.en_ci = None  # excitation energies from external code
        self.Vijn = np.array([None])  # potential from external code if PCM
        self.Vijn_fortran_flip =  None
        ### following parameters are computed by means of an external quantum chemistry code, if needed ####
        self.hamiltonian = None
        self.hamiltonian_t = None
        self.ao_coefficients = None
        self.two_electron_int = None
        self.ao_kinetic_integral = None
        self.ao_potential_integral = None
        self.doubly_occupied_mo = None  #### at the moment we consider only simple cases of singlet states
        self.num_mo = None
        self.control_operator = None ## questa terminologia non è molto appropriata per una molecola 
                                     ## (qui ci riferiamo nel caso specifico all'operatore: \sum_{p,q} <p|\sum_i Z_i/r_i|q>)
                                     ## sulla base dei determinanti di Slater
        self.basis_determinants = None