import numpy as np

from read_and_set.read import auxiliary_functions as af

from propagator.ABCPropagator import ABCPropagator
from propagator.PropagatorTerms import PropagatorTerms

from SystemObj import Func_tMatrix

#import math_functions as mf

class PropagatorOCfwd(ABCPropagator):
    def __init__(self):
        super().__init__()
        self.propagator_terms = PropagatorTerms()
        self.propagator = []

    def set_propagator(self, molecule, env):
        self.init_propagaror_terms(molecule, env)
        self.clean_propagator()
        self.add_term_to_propagator("eulero1_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if self.propagator_terms.pcm != None:
            if self.propagator_terms.pcm.par.env == "sol":
                self.add_term_to_propagator("eulero_pcm")


    def propagate_one_step(self, i, dt, field_dt_vector):
        for func in self.propagator:
            func(i, 1, dt, field_dt_vector)


    def propagate_n_step(self, discrete_time_par, field):
        if((field.time_axis[1] - field.time_axis[0]) - discrete_time_par.dt > discrete_time_par.dt * 0.001):
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



class PropagatorOCbwd(ABCPropagator):
    def __init__(self):
        super().__init__()
        self.propagator_terms = PropagatorTerms()
        self.propagator = []

    def set_propagator(self, molecule, env):
        self.init_propagaror_terms(molecule, env)
        self.clean_propagator()
        self.add_term_to_propagator("eulero1_coeff")
        self.add_term_to_propagator("eulero_energy")
        self.add_term_to_propagator("eulero_field")
        if self.propagator_terms.pcm != None:
            if (self.propagator_terms.pcm.par.env == "sol"):
                self.add_term_to_propagator("oc_pcm_bwd")


    def propagate_one_step(self, i, dt, field_dt_vector, wf_fwd):
        for func in self.propagator:
            func(i, 1, -dt, field_dt_vector, wf_fwd)

    def propagate_n_step(self, discrete_time_par, field, wf_fwd):
        if((field.time_axis[1] - field.time_axis[0]) - discrete_time_par.dt > discrete_time_par.dt *0.001):
            af.exit_error("ERROR: Propagation time step and field time step are different")
        wf_matrix_out = Func_tMatrix()
        wf_matrix_out.time_axis = np.linspace(0,
                                              discrete_time_par.dt * discrete_time_par.nstep,
                                              discrete_time_par.nstep + 1)
        out = list()
        out.append(self.propagator_terms.mol.wf.ci)
        for i in range(discrete_time_par.nstep):
            self.propagate_one_step(i, discrete_time_par.dt, field.f_xyz[i], wf_fwd)
            out.append(self.propagator_terms.mol.wf.ci)
        wf_matrix_out.f_xyz = np.array(out)
        return wf_matrix_out










