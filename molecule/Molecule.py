from molecule.WaveFunction import WaveFunction
from parameters.MoleculeParameters import MoleculeParameters
from read_and_set.read import auxiliary_functions as af
from molecule.Psi4_DetHam_constructor import Psi4InterfaceMolecule


class Molecule:
    def __init__(self):
        self.wf = WaveFunction()
        self.par = MoleculeParameters()

    def init_molecule(self, molecule_input):
        self.wf.set_wf(molecule_input.wf_ci, True)
        self.par.muT = molecule_input.muT
        self.par.en_ci = molecule_input.en_ci
        self.par.Vijn = molecule_input.Vijn
        if(self.par.Vijn.all() != None):
            self.par.Vijn_fortran_flip = af.flip_3D_py2f(self.par.Vijn)
        self.par.num_mo = molecule_input.nmo
        self.par.doubly_occupied_mo = molecule_input.doubly_occupied_mo
        self.par.ao_coefficients = molecule_input.ao_coefficients
        self.par.ao_kinetic_integral = molecule_input.ao_kinetic_integral
        self.par.ao_potential_integral = molecule_input.ao_potential_integral
        self.par.two_electron_int = molecule_input.two_electron_int
        
    def init_molecular_hamiltonian(self):
        self.set_molecular_hamiltonian_configuration_basis()
        self.set_control_operator_electron_nuclei()
        
    def set_n_electrons_basis_set(self):
        self.par.basis_determinants = Psi4InterfaceMolecule.set_n_electrons_basis_set(self.par)
        
    def set_molecular_hamiltonian_configuration_basis(self):
        self.set_n_electrons_basis_set()
        self.par.hamiltonian = Psi4InterfaceMolecule.molecular_hamiltonian_configuration_basis(self.par)
        
    def set_control_operator_electron_nuclei(self):
        self.par.control_operator = Psi4InterfaceMolecule.electron_nuclei_determinant_basis(self.par)
        
    def propagate_hamiltonian(self, perturbation):
        self.par.hamiltonian_t = self.par.hamiltonian - (1-perturbation)*self.par.control_operator





























































