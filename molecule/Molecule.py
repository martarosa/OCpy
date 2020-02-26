from molecule.WaveFunction import WaveFunction
from parameters.MoleculeParameters import MoleculeParameters
from read_and_set.read import auxiliary_functions as af


class Molecule:
    def __init__(self):
        self.wf = WaveFunction()
        self.par = MoleculeParameters()

    def init_molecule(self, molecule_input):
        self.wf.set_wf(molecule_input.wf_ci, True)
        self.par.muT = molecule_input.muT
        self.par.en_ci = molecule_input.en_ci
        self.par.Vijn = molecule_input.Vijn
        if(self.par.Vijn != None):
            self.par.Vijn_fortran_flip = af.flip_3D_py2f(self.par.Vijn)



























































