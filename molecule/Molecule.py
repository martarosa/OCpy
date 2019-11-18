from molecule.WaveFunction import WaveFunction
from read import auxiliary_functions as af

class Molecule:
    def __init__(self):
        self.wf = WaveFunction()
        self.muT = None  # transition dipoles from external code
        self.en_ci = None  # excitation energies from external code
        self.Vijn = None   # potential from external code if PCM
        self.Vijn_fortran_flip = None


    def init_molecule(self, molecule_parameters):
        self.wf.set_wf(molecule_parameters.wf_ci, True)
        self.muT = molecule_parameters.muT
        self.en_ci = molecule_parameters.en_ci
        self.Vijn = molecule_parameters.Vijn
        if(self.Vijn.any()):
            self.Vijn_fortran_flip = af.flip_3D_py2f(self.Vijn)







