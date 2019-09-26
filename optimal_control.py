###da sistemare rabitz uno e due e che legge omega anche quando ha i campi dove non serve
##### fare la versione fortran aggiornata a qt

import numpy as np
import sys

import PropagatorOCRabitz as rabitzI
import PropagatorsEulero as prop
import auxiliary_functions as af

from Propagator import Propagator
from Field import PropagatorFieldOC, Field
from SaveOC import SaveOC
#from genetic import GeneticParameters, Genetic, Recombination, Mutation

class OC:
    def __init__(self):
        self.nstep = None
        self.dt = None
        self.oc_algorithm_name = None
        self.oc_iterator = None

        #self.prop_psi = Propagator()
        self.psi_coeff_t = None
        self.field_vector_t = None

        self.state_to_project_on = None
        self.initial_c0 = None

        self.save = SaveOC()
        self.parameters = dict({})
        self.alpha_t = None
        self.current_iteration = 0
        self.J = None


    def initialize_oc(self, input_par_file, mol, field, env):
        self.initialize_common_oc_parameters(input_par_file)
        self.oc_iterator.init_oc(mol, field, env)
        self.save = self.oc_iterator.init_save(self.save, input_par_file)



    def initialize_propagation_attributes(self, molecule, field, env):
        pass

    def iterate_oc(self):
        pass


    def calc_every_iteration_output(self):
        pass

    def set_initial_c0(self):
        pass

    def calc_J(self):
        pass


    def initialize_common_oc_parameters(self, input_par_file):
        self.nstep = int(input_par_file.get_sys_par()["nstep"])
        self.dt = float(input_par_file.get_sys_par()["dt"])
        #self.restart = input_par_file.get_oc_par()["restart"]
        #self.oc_algorithm = input_par_file.get_sys_par()["propagation"]  # questo decide il propagatore genetico, rabitz, eulero ecc


        parameters = ['alpha', 'alpha0', 'convergence_thr','n_iterations', 'Ns', 'delta_ts']
        values = [ input_par_file.get_oc_par()["alpha"],
                   float(input_par_file.get_oc_par()["alpha0"]),
                   float(input_par_file.get_oc_par()["convergence_thr"]),
                   int(input_par_file.get_oc_par()["n_iterations"]),
                   float(input_par_file.get_oc_par()["ns"]),
                   float(input_par_file.get_oc_par()["delta_ts"])]
        self.parameters = dict(zip(parameters, values))


        self.alpha_t = np.zeros(self.nstep)
        self.set_alpha_t()
        self.J = 99999
        self.convergence_t = 99999
        self.state_to_project_on = af.normalize_vector(
            [float(i) for i in input_par_file.get_oc_par()["state_to_project_on"].split(' ')])

        self.save.set_save(input_par_file.get_sys_par()["folder"], input_par_file.get_sys_par()["name"])
        self.save.set_n_save([int(input_par_file.get_save_par()["save_iteration_pop_step"]),
                              int(input_par_file.get_save_par()["save_iteration_field_step"])],
                             int(input_par_file.get_save_par()["save_restart_step"]))




    def set_alpha_t(self):
        for i in range(self.nstep):
            if self.parameters['alpha'] == "const":
                self.alpha_t[i] = self.parameters['alpha0']  # const
            elif self.parameters['alpha'] == "sin":
                self.alpha_t[i] = self.parameters['alpha0'] / np.square(np.sin(np.pi * (i + 1) / self.nstep))  # paper gross
            elif self.parameters['alpha'] == "quin":
                self.alpha_t[i] = self.parameters['alpha0'] / np.exp(
                    -np.power((((i) * self.dt - 125) / 220), 12))  # paper quinolone


    def set_starting_field(self, starting_field):
        pass

    def check_convergence(self):
        J_prev_tmp = np.copy(self.J)
        self.calc_J()
        self.convergence_t = self.J - J_prev_tmp

    def get_dt(self):
        return self.dt


    def calc_proj(self):
        projection_on_target_state = af.projector_mean_value(self.prop_psi.mol.wf.Ci(), self.state_to_project_on)\
                                     /(np.dot(self.prop_psi.mol.wf.Ci(), np.conj(self.prop_psi.mol.wf.Ci())))
        return np.real(projection_on_target_state)


    def calc_field_integral(self):
        integral_field = af.field_J_integral(self.field_vector_t, self.dt)
        return np.real(integral_field)





class OC_Genetic(OC):
    def __init__(self):
        super().__init__()
        self.genetic = GeneticParameters()


    def initialize_oc(self, input_par_file, mol, field, env):
        self.initialize_common_oc_parameters(input_par_file)
        self.initialize_propagation_attributes(mol, field, env)
        #pop_t e field_t tengono il risultato migliore
        #field_par tiene (tutte? le n per ripartire?) repliche dei parametri
        #log tiene la convergenza, la fitness,
        self.save.creation_save_files(
                                 [".log", "_final_pop.dat", "_pop_t.dat", "_field_t.dat"],
                                 ["#Final states populations \n #fields: n_iteration, convergence, J, projection on target state, alpha*integral(field^2)\n" ,
                                  "#Final states populations \n #fields: n_iteration, states population  \n",
                                  "#Iteration populations \n #fields: n_iteration, time, states population \n",
                                  "#Optimized field \n#fields: n_iteration, time, fields components  \n"],
                                 input_par_file)



    def initialize_propagation_attributes(self, molecule, field, env):
        prop_psi = Propagator()
        prop_psi.set_propagator(self.dt, molecule, env)
        self.prop_psi = np.ndarray(shape=(self.genetic.n_genes), dtype=Propagator())
        for i in range(self.genetic.n_genes):
            self.prop_psi[i] = prop_psi
        self.psi_coeff_t = np.zeros([self.genetic.n_genes, self.nstep + 1, self.prop_psi[0].mol.wf.n_ci], dtype=complex)
        self.J = np.zeros([self.genetic.n_genes], dtype=complex)
        self.set_initial_c0()


    def iterate_oc(self):
        for i in range(self.genetic.n_genes):
            self.prop_psi[i].propagate_n_step(self.nstep, self.prop_field[i])

        while (self.current_iteration <= self.parameters["n_iterations"] or self.parameters["convergence_thr"] < self.convergence_t):
            self.check_convergence()
            #do genetic
        for i in range(self.genetic.n_genes):
            self.prop_psi[i].mol.wf.set_ci(self.initial_c0, 0)
            self.psi_coeff_t[i] = self.prop_psi[i].propagate_n_step(self.nstep, self.prop_field[i])
            self.save.save_every_iteration(self.current_iteration, [".log"],
                                   np.real(self.calc_every_iteration_output()))


        #self.save.save_every_n_iterations(self.current_iteration, ["_pop_t.dat", "_field_t.dat"],
        #                              [np.real(af.population_from_wf_matrix(self.psi_coeff_t)), self.field_vector_t],
        #                              "_field_bkp.dat",
        #                              self.field_vector_t)



    def calc_J(self):
        for i in range(self.genetic.n_genes):
            self.J[i] = np.real(af.projector_mean_value(self.prop_psi[i].mol.wf.Ci(), self.state_to_project_on) \
                 - af.alpha_field_J_integral(self.field_vector_t[i], self.alpha_t, self.dt))


    def set_starting_field(self, starting_field):
        self.field_vector_t = starting_field


    def calc_every_iteration_output(self):
        out = np.array([self.convergence_t[:KEEP], self.J[:KEEP], self.calc_proj()[:KEEP], self.calc_field_integral()[:KEEP]])
        #pop = np.zeros([self.genetic.n_genes])
        #for i in range(self.genetic.n_genes):
        #    pop[i] = af.population_from_wf_vector(self.prop_psi[i].mol.wf.get_ci())
        #return np.array([out.T, np.real(pop)])
        return out.T


    def set_initial_c0(self):
        self.initial_c0 = self.prop_psi[0].mol.wf.Ci()
        if self.state_to_project_on.size != self.initial_c0.size:
            sys.exit(
                "ERROR: Wrong dimension of target state: number of coefficients is different from wavefunction. ")



    def set_genetic_fields(self, field_type, t0):
        field_template = Field()
        self.prop_field = np.ndarray(shape=(self.genetic.n_genes), dtype=Field)
        for i in range(self.genetic.n_genes):
            self.prop_field[i] = field_template
            self.prop_field[i].set_field_parameters(field_type, self.genetic.ampitudes_cromosomes, 0, 0, self.genetic.omega, t0)
            self.prop_field[i].create_field(self.nstep, self.dt)