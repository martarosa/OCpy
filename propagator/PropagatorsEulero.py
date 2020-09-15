import numpy as np

from medium.ABCMedium import ABCMedium
from molecule.Molecule import Molecule

from read_and_set.read import auxiliary_functions as af

from propagator.ABCPropagator import ABCPropagator
from propagator.ClassicalPropagatorTerms import ClassicalPropagatorTerms

from SystemObj import Func_tMatrix



class PropagatorEulero1Order(ABCPropagator):
    def __init__(self):
        super().__init__()
        self.mol = Molecule()
        self.medium = ABCMedium()
        self.propagator_terms = ClassicalPropagatorTerms()
        self.propagator = []


    def set_propagator(self, molecule, medium):
        self.mol = molecule
        self.medium = medium
        self.clean_propagator()
        self.propagator_terms.init()
        self.add_term_to_propagator("eulero1_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if self.medium != None:
            if self.medium.par.medium == "sol":
                self.add_term_to_propagator("eulero_medium")
        self.add_term_to_propagator("norm")


    def propagate_one_step(self, i, dt, field_dt_vector):
        for func in self.propagator:
            func(i, 1, dt, self.mol, field_dt_vector, self.medium)


    def propagate_n_step(self, discrete_time_par, field):
        if((field.time_axis[1] - field.time_axis[0]) - discrete_time_par.dt > discrete_time_par.dt *0.001):
            af.exit_error("ERROR: Propagation time step and field time step are different")
        wf_matrix_out = Func_tMatrix()
        wf_matrix_out.time_axis = np.linspace(0,
                                              discrete_time_par.dt * discrete_time_par.nstep,
                                              discrete_time_par.nstep + 1)
        out = list()
        out.append(self.mol.wf.ci)
        for i in range(discrete_time_par.nstep):
            self.propagate_one_step(i, discrete_time_par.dt, field.f_xyz[i])
            out.append(self.mol.wf.ci)
        wf_matrix_out.f_xyz = np.array(out)
        return wf_matrix_out




class PropagatorEulero2Order(ABCPropagator):
    def __init__(self):
        super().__init__()
        self.mol = Molecule()
        self.medium = ABCMedium()
        self.propagator_terms = ClassicalPropagatorTerms()
        self.propagator = []

    def set_propagator(self, molecule, medium):
        self.mol = molecule
        self.medium = medium
        self.clean_propagator()
        self.propagator_terms.init()
        self.add_term_to_propagator("eulero2_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if self.medium != None:
                self.add_term_to_propagator("eulero_medium")
        self.add_term_to_propagator("norm")


    def propagate_one_step(self, dt, field_dt_vector, order=2):
        for func in self.propagator:
            func(order, dt, self.mol, field_dt_vector, self.medium)

    def propagate_n_step(self, discrete_time_par, field):
        if((field.time_axis[1] - field.time_axis[0]) - discrete_time_par.dt > discrete_time_par.dt *0.001):
            af.exit_error("ERROR: Propagation time step and field time step are different")
        wf_matrix_out = Func_tMatrix()
        wf_matrix_out.time_axis = np.linspace(0,
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
        wf_matrix_out.f_xyz = np.array(out)
        return wf_matrix_out