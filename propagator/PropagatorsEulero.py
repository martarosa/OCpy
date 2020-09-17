import numpy as np

from read_and_set.read import auxiliary_functions as af

from propagator.ABCPropagator import ABCPropagator
from propagator.PropagatorTerms import PropagatorTerms

from SystemObj import Func_tMatrix



class PropagatorEulero1Order(ABCPropagator):
    def __init__(self):
        super().__init__()
        self.propagator_terms = PropagatorTerms()
        self.propagator = []


    def set_propagator(self, molecule, medium):
        self.init_propagator_terms(molecule, medium)
        self.clean_propagator()
        self.add_term_to_propagator("eulero1_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if self.propagator_terms.medium != None:
            if self.propagator_terms.medium.par.medium == "sol":
                self.add_term_to_propagator("eulero_medium")
        self.add_term_to_propagator("norm")


    def propagate_one_step(self, i, dt, field_dt_vector):
        for func in self.propagator:
            func(i, 1, dt, field_dt_vector)


    def propagate_n_step(self, discrete_time_par, field):
        if((field.time_axis[1] - field.time_axis[0]) - discrete_time_par.dt > discrete_time_par.dt *0.001):
            af.exit_error("ERROR: Propagation time step and field time step are different")
        wf_matrix_out = Func_tMatrix()
        wf_matrix_out.time_axis = np.linspace(0,
                                              discrete_time_par.dt * discrete_time_par.nstep,
                                              discrete_time_par.nstep + 1)
        out = list()
        out.append(self.propagator_terms.mol.wf.ci)
        for i in range(discrete_time_par.nstep):
            self.propagate_one_step(i, discrete_time_par.dt, field.f_xyz[i])
            out.append(self.propagator_terms.mol.wf.ci)
        wf_matrix_out.f_xyz = np.array(out)
        return wf_matrix_out






class PropagatorEulero2Order(ABCPropagator):
    def __init__(self):
        super().__init__()
        self.propagator_terms = PropagatorTerms()
        self.propagator = []

    def set_propagator(self, molecule, medium):
        self.init_propagator_terms(molecule, medium)
        self.clean_propagator()
        self.add_term_to_propagator("eulero2_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if self.propagator_terms.medium != None:
            #if self.propagator_terms.medium.par.medium == "sol":
                self.add_term_to_propagator("eulero_medium")
        self.add_term_to_propagator("norm")


    def propagate_one_step(self, dt, field_dt_vector, order=2):
        for func in self.propagator:
            func(order, dt, field_dt_vector)

    def propagate_n_step(self, discrete_time_par, field):
        if self.propagator_terms.medium.par.medium == "nanop":
            self.propagator_terms.medium.reset_medium(self.propagator_terms.mol, field)
        if((field.time_axis[1] - field.time_axis[0]) - discrete_time_par.dt > discrete_time_par.dt *0.001):
            af.exit_error("ERROR: Propagation time step and field time step are different")
        wf_matrix_out = Func_tMatrix()
        wf_matrix_out.time_axis = np.linspace(0,
                                              discrete_time_par.dt * discrete_time_par.nstep,
                                              discrete_time_par.nstep + 1)
        out = list()
        out.append(self.propagator_terms.mol.wf.ci)
        for i in range(discrete_time_par.nstep):
            #print("python field")
            #print(i+1)
            #print(field.f_xyz[i][2])
            if i != 0:
                self.propagate_one_step(discrete_time_par.dt, field.f_xyz[i])
            else:
                self.propagate_one_step(discrete_time_par.dt, field.f_xyz[i], order=1)
            out.append(self.propagator_terms.mol.wf.ci)
        wf_matrix_out.f_xyz = np.array(out)
        return wf_matrix_out