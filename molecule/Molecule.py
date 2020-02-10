from molecule.WaveFunction import WaveFunction

class MoleculeParameters:
    def __init__(self):
        self.muT = None  # transition dipoles from external code
        self.en_ci = None  # excitation energies from external code
        self.Vijn = None  # potential from external code if PCM


class Molecule:
    def __init__(self):
        self.wf = WaveFunction()
        self.par = MoleculeParameters()

    def init_molecule(self, molecule_input):
        self.wf.set_wf(molecule_input.wf_ci, True)
        self.par.muT = molecule_input.muT
        self.par.en_ci = molecule_input.en_ci
        self.par.Vijn = molecule_input.Vijn



























































