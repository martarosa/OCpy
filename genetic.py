import numpy as np
import sys

import PropagatorOCRabitz as rabitzI
import PropagatorsEulero as prop
import auxiliary_functions as af

from Propagator import Propagator
from Field import PropagatorFieldOC, Field
from environ import PCM
from SaveOC import SaveOC

import numpy as np



class GeneticParameters:
    def __init__(self):
        self.omega = None
        self.n_omega = None
        #self.amplitudes_gene = None
        self.n_genes = None
        self.ampitudes_cromosomes = None
        self.n_cromosomes = None


    def set_parameters(self, delta, n_cromosomi, mol):
        self.n_omega = mol.en_ci[-1]/delta
        self.n_genes = self.n_omega
        self.n_cromosomes = n_cromosomi
        self.omega = np.zeros(self.n_omega)
        self.omega[0] = delta
        for i in range(self.n_omega):
            self.omega[i + 1] = self.omega[i] + delta


    def init_first_generation_amplitudes(self, low, high):
        self.ampitudes_cromosomes = np.random.uniform(low=low, high=high, size=(self.n_genes, self.n_cromosomes))




class Genetic():
    def __init__(self):



class Recombination():
    def __init__(self):
        super().__init__()


    def evolve(self):



class Mutation():
    def __init__(self):
        super().__init__()







class OptimalControl_genetic:
    def __init__(self):
        self.prop_psi = None
        self.fields = None
        self.genetic = GeneticParameters()


        self.save = SaveOC()
        self.psi_coeff_t = None

        self.field_vector_t = None
        self.nstep = None
        self.dt = None
        self.restart = None
        self.convergence = None
        self.convergence_thr = None
        self.convergence_t = None
        self.state_to_project_on = None
        self.n_iterations = None
        self.initial_c0 = None
        self.current_iteration=0
        #output parameters
        self.folder = None                      #input and output folder
        self.name_output = None                 #prefix of output files
        self.J = None                           #functional to minimize in optimal control
        self.J_prev = None
        self.alpha = None
        self.alpha0 = None
        self.alpha_t = None


    def set_fields(self, field_type, t0):
        field_template = Field()
        self.fields = np.ndarray(shape=(self.genetic.n_genes), dtype=Field)
        for i in range(self.genetic.n_genes):
            self.fields[i] = field_template
            self.fields[i].set_field_parameters(field_type, self.genetic.ampitudes_cromosomes, 0, 0, self.genetic.omega, t0)
            self.fields[i].create_field(self.nstep, self.dt)


    def initialize_propagation(self, input_par_file, dt, molecule, env):
        prop_psi = Propagator()
        prop_psi.set_propagator(dt, molecule, env)
        if self.state_to_project_on.size != self.initial_c0.size:
            sys.exit("ERROR: Wrong dimension of target state: number of coefficients is different from wavefunction. ")
        self.prop_psi = np.ndarray(shape=(self.genetic.n_genes), dtype=Propagator())
        for i in range(self.genetic.n_genes):
            self.prop_psi[i] = prop_psi


    def iterate_oc(self):
        for i in range(self.genetic.n_genes):
            self.prop_psi[i].propagate_n_step(self.nstep, self.fields[i])

        while (self.current_iteration <= self.n_iterations or self.convergence_thr < self.convergence_t):
            self.check_convergence()
            #do genetic
        for i in range(self.genetic.n_genes):
            self.prop_psi[i].propagate_n_step(self.nstep, self.fields[i])










    def check_convergence(self):
        self.J_prev = np.copy(self.J)
        for i in range(self.genetic.n_genes):
            self.J[i] = af.projector_mean_value(self.prop_psi[i].mol.wf.Ci(),
                                                self.state_to_project_on)\
                        -af.alpha_field_J_integral(self.field[i],
                                          self.alpha_t,
                                          self.dt)
        #print(self.current_iteration, self.J)
        #self.convergence_t=self.J-self.J_prev





