import numpy as np

import PropagatorsEulero as prop
import auxiliary_functions as af


from OCIterator import OCIterator
from Field import Field

class GeneticParameters:
    def __init__(self):
        self.n_amplitudes_omega_genes = None #n_ampitudes = n_omega
        self.delta_omega = None
        self.genes_amplitude_limits = []



class InitGeneticPar():
    def __init__(self):
        self.genetic_parameters = GeneticParameters()

    def init(self, user_input):
        self.genetic_parameters.n_amplitudes_omega_genes = user_input.chr.par['n_genes']
        self.genetic_parameters.genes_amplitude_limits.append(user_input.chr.par['amplitude_low_limit'])
        self.genetic_parameters.genes_amplitude_limits.append(user_input.chr.par['amplitude_hight_limit'])
        self.genetic_parameters.delta_omega = user_input.chr.par['delta_omega']




class Chromosome():
    def __init__(self):
        self.n_amplitude_genes = None #n_ampitudes
        self.amplitude_genes_vector = None #amplitudes
        self.field = Field()
        self.psop_psi = prop.PropagatorEulero2Order()


    def init_chromosome(self, chromosome_parameters,  molecule, starting_field, pcm, oc_parameters, alpha_t):
        self.n_amplitude_genes = chromosome_parameters.n_genes

        self.set_omega_vector(molecule)










class OCGeneticIterator(OCIterator):
    def __init__(self):
        super().__init__()
        self.nstep = None
        self.dt = None
        self.target_state = None
        self.alpha_t = None
        self.convergence_t =  99999
        self.J = 99999

        self.field_psi_matrix = None
        self.psi_coeff_t = None

        self.dict_out = {}
        self.dict_restart = {}

        #genetic

        self.n_chromosomes = None
        self.n_amplitudes_omega_genes = None
        self.genes_amplitude_limits = []
        self.omegas_vector = None
        self.delta = None



    def iterate(self, current_iteration):
        self.psi_coeff_t = self.prop_psi.propagate_n_step(self.nstep, self.field_psi_matrix)


    def init_oc_iterator(self, molecule, starting_field, env, oc_parameters, alpha_t):
        self.nstep = oc_parameters.nstep
        self.dt = oc_parameters.dt
        self.target_state = oc_parameters.target_state
        self.alpha_t = alpha_t
        self.n_amplitudes_omega_genes = oc_pa



        self.prop_psi.set_propagator(self.dt, molecule, env)
        self.field_psi_matrix = np.copy(starting_field.field)


    def set_omega_vector(self, molecule):
        self.delta = molecule.en_ci[-1]/self.n_amplitude_genes
        self.omegas_vector = np.zeros(self.n_amplitude_genes)
        self.omegas_vector[0] = self.delta
        for i in range(self.n_amplitude_genes):
            self.omegas_vector[i + 1] = self.omega_vector[i] + self.delta




        self.init_output_dictionary()





    def init_output_dictionary(self):
        self.dict_out['pop_t'] = self.get_pop_t


    def get_pop_t(self):
        pop_t = np.real(af.population_from_wf_matrix(self.psi_coeff_t))
        return pop_t

    def get_restart(self):
        return self.field_psi_matrix