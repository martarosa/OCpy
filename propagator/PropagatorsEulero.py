import numpy as np


from read_and_set.read import auxiliary_functions as af

from propagator.ABCPropagator import ABCPropagator
from molecule.Molecule import Molecule
from SystemObj import Func_tMatrix
import dictionaries.PropagatorTermsDictionaries as ptdict




class PropagatorEulero1Order(ABCPropagator):
    def __init__(self):
        super().__init__()
        self.mol = Molecule()
        self.medium = None
        self.propagator_terms = None
        self.propagator = []
        self.wf_matrix_out = Func_tMatrix()

    def init(self, molecule, medium, propagator):
        self.mol = molecule
        self.medium = medium
        self.propagator_terms = ptdict.PropagatorTermsDict[propagator]()
        self.propagator_terms.init()
        
    def clean_propagator(self):
        self.propagator = []


    def set_propagator(self, molecule, medium):
        self.init(molecule, medium, "eulero_1order")
        self.clean_propagator()
        self.add_term_to_propagator("eulero1_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if self.medium != None:
            if self.medium.par.medium == "sol":
                self.add_term_to_propagator("eulero_medium")
        self.add_term_to_propagator("norm")


    def propagate_one_step(self, dt, field_dt_vector):
        for func in self.propagator:
            func(self.mol, 1, dt, field_dt_vector, self.medium)


    def propagate_n_step(self, discrete_time_par, field):
        if((field.time_axis[1] - field.time_axis[0]) - discrete_time_par.dt > discrete_time_par.dt *0.001):
            af.exit_error("ERROR: Propagation time step and field time step are different")

        self.wf_matrix_out.time_axis = np.linspace(0,
                                              discrete_time_par.dt * discrete_time_par.nstep,
                                              discrete_time_par.nstep + 1)
        out = list()
        out.append(self.mol.wf.ci)
        for i in range(discrete_time_par.nstep):
            self.propagate_one_step(discrete_time_par.dt, field.f_xyz[i])
            out.append(self.mol.wf.ci)
        self.wf_matrix_out.f_xyz = np.array(out)





class PropagatorEulero2Order(ABCPropagator):
    def __init__(self):
        super().__init__()
        self.mol = Molecule()
        self.medium = None

        self.propagator_terms = None
        self.propagator = []

    def init(self, molecule, medium, propagator):
        self.mol = molecule
        self.medium = medium
        self.propagator_terms = ptdict.PropagatorTermsDict[propagator]()
        self.propagator_terms.init()
        self.wf_matrix_out = Func_tMatrix()

    def clean_propagator(self):
        self.propagator = []

    def set_propagator(self, molecule, medium, prop_conf = None):
        self.init(molecule, medium, "eulero_2order")
        self.clean_propagator()
        self.add_term_to_propagator("eulero2_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if self.medium != None:
            self.add_term_to_propagator("eulero_medium")
        self.add_term_to_propagator("norm")


    def propagate_one_step(self, dt, field_dt_vector, order=2):
        for func in self.propagator:
            func(self.mol, order, dt, field_dt_vector, self.medium)

    def propagate_n_step(self, discrete_time_par, field):
        if(field.time_axis[1] - field.time_axis[0]) - discrete_time_par.dt > discrete_time_par.dt *0.001:
            af.exit_error("ERROR: Propagation time step and field time step are different")

        self.wf_matrix_out.time_axis = np.linspace(0,
                                              discrete_time_par.dt * discrete_time_par.nstep,
                                              discrete_time_par.nstep + 1)
        out = list()
        out.append(self.mol.wf.ci)
        for i in range(discrete_time_par.nstep):

            if i != 0:
                self.propagate_one_step(discrete_time_par.dt, field.f_xyz[i])
            else:
                self.propagate_one_step(discrete_time_par.dt, field.f_xyz[i], order=1)
            out.append(self.mol.wf.ci)
        self.wf_matrix_out.f_xyz = np.array(out)
        
class PropagatorEulero2OrderGeneralPerturbation(ABCPropagator):
    def __init__(self):
        super().__init__()
        self.mol = Molecule()
        self.medium = None
        
        self.propagator_terms = None
        self.propagator = []
        
    def init(self, molecule, medium, propagator):
        self.mol = molecule
        self.medium = medium
        self.propagator_terms = ptdict.PropagatorTermsDict[propagator]()
        self.propagator_terms.init()
        self.wf_matrix_out = Func_tMatrix()
        
    def clean_propagator(self):
        self.propagator = []
        
    def set_propagator(self, molecule, medium, prop_conf = None):
        self.init(molecule, medium, "eulero_2order_psi4")
        self.clean_propagator()
        self.add_term_to_propagator("eulero2_coeff")
        if af.is_diagonal(self.mol.par.hamiltonian):
            self.add_term_to_propagator("eulero_energy")
        else:
            self.add_term_to_propagator("eulero_external_hamiltonian")
      #  self.add_term_to_propagator("eulero_perturbation")
        if self.medium != None:
            self.add_term_to_propagator("eulero_medium")
        self.add_term_to_propagator("norm")
    
    def propagate_one_step(self, dt, perturbation, order=2):
        for func in self.propagator:
            func(self.mol, order, dt, perturbation, self.medium)
            
    def propagate_n_step(self, discrete_time_par, field):
        if(field.time_axis[1] - field.time_axis[0]) - discrete_time_par.dt > discrete_time_par.dt *0.001:
            af.exit_error("ERROR: Propagation time step and field time step are different")

        self.wf_matrix_out.time_axis = np.linspace(0,
                                              discrete_time_par.dt * discrete_time_par.nstep,
                                              discrete_time_par.nstep + 1)
        out = list()
        out.append(self.mol.wf.ci)
        for i in range(discrete_time_par.nstep):

            if i != 0:
                self.propagate_one_step(discrete_time_par.dt, field.f_xyz[i])
            else:
                self.propagate_one_step(discrete_time_par.dt, field.f_xyz[i], order=1)
            out.append(self.mol.wf.ci)
        self.wf_matrix_out.f_xyz = np.array(out)