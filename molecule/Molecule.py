from molecule.WaveFunction import WaveFunction
from read import auxiliary_functions as af

class Molecule:
    def __init__(self):
        self.wf = WaveFunction()
        self.muT = None  # transition dipoles from external code
        self.en_ci = None  # excitation energies from external code
        self.Vijn = None   # potential from external code if PCM
        self.Vijn_fortran_flip = None


    def init_molecule(self, molecule_input):
        self.wf.set_wf(molecule_input.wf_ci, True)
        self.muT = molecule_input.muT
        self.en_ci = molecule_input.en_ci
        self.Vijn = molecule_input.Vijn
        if(self.Vijn!= None):
            self.Vijn_fortran_flip = af.flip_3D_py2f(self.Vijn)

























































