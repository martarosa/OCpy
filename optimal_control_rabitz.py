import numpy as np
import sys

import PropagatorOCRabitz as rabitzI
import PropagatorsEulero as prop
import auxiliary_functions as af

from Propagator import Propagator
from Field import PropagatorFieldOC, Field
from SaveOC import SaveOC
from optimal_control import OC


class OC_Rabitz(OC):
    def __init__(self):
        super().__init__()

        self.prop_psi = rabitzI.PropagatorOCfwd()
        self.prop_chi = rabitzI.PropagatorOCbwd()
        self.prop_field = PropagatorFieldOC()

        # storage vector optimal control
        self.chi_coeff_t = None
        self.field_chi_vector_t = None
        self.convergence = None
        self.convergence_t = None
        self.oc_algorithm = "rabitzi"





    def initialize_oc(self, input_par_file, mol, field, env):
        self.initialize_common_oc_parameters(input_par_file)
        self.initialize_propagation_attributes(mol, field, env)
        print("quiqui")
        print(input_par_file.sys.pa)
        self.save.creation_save_files(
                                 [".log", "_final_pop.dat", "_pop_t.dat", "_field_t.dat"],
                                 ["#Final states populations \n #fields: n_iteration, convergence, J, projection on target state, alpha*integral(field^2)\n" ,
                                  "#Final states populations \n #fields: n_iteration, states population  \n",
                                  "#Iteration populations \n #fields: n_iteration, time, states population \n",
                                  "#Optimized field \n#fields: n_iteration, time, fields components  \n"],
                                 input_par_file)


    def initialize_propagation_attributes(self, molecule, field, env):
        self.prop_psi.set_propagator(self.dt, molecule, env)
        self.prop_chi.set_propagator(-1 * self.dt, molecule, env)
        self.set_starting_field(field)
        self.psi_coeff_t = np.zeros([self.nstep + 1, self.prop_psi.mol.wf.n_ci], dtype=complex)
        self.psi_coeff_t[0] =  self.prop_psi.mol.wf.Ci()
        self.chi_coeff_t = np.zeros([self.nstep + 1, self.prop_psi.mol.wf.n_ci], dtype=complex)
        self.field_chi_vector_t = np.zeros([self.nstep, 3], dtype=complex)
        self.set_initial_c0()




    def calc_every_iteration_output(self):
        out = [self.convergence_t, self.J, self.calc_proj(), self.calc_field_integral()]
        pop = af.population_from_wf_vector(self.prop_psi.mol.wf.Ci())
        return [out, np.real(pop)]




    def calc_J(self):
        if self.oc_algorithm == "rabitzi":
            self.J = np.real(af.projector_mean_value(self.prop_psi.mol.wf.Ci(), self.state_to_project_on) \
                     - af.alpha_field_J_integral(self.field_vector_t, self.alpha_t, self.dt))
        elif self.oc_algorithm == "rabitzii":
            self.J = np.real(2 * np.real(np.dot(self.state_to_project_on, self.prop_psi.mol.wf.Ci()) \
                                 - af.alpha_field_J_integral(self.field_vector_t, self.alpha_t, self.dt)))



    def iterate_oc(self):
        self.save.save_every_iteration(-1, [".log", "_final_pop.dat"], np.real(self.calc_every_iteration_output()))
        while (self.current_iteration <= self.parameters["n_iterations"] or self.parameters["convergence_thr"] < self.convergence_t):
            if self.current_iteration != 0:
                for i in range(self.nstep, 0, -1):
                    if i == self.nstep:
                        if self.oc_algorithm == 'rabitzi':
                            self.prop_chi.mol.wf.set_ci(af.apply_projection(self.prop_psi.mol.wf.Ci(),
                                                                            self.state_to_project_on),
                                                        0)
                        elif self.oc_algorithm == 'rabitzii':
                            self.prop_chi.mol.wf.set_ci(self.state_to_project_on, 0)


                        self.chi_coeff_t[i] = self.prop_chi.mol.wf.Ci()
                    self.field_chi_vector_t[i - 1] = self.prop_field.propagate_field_OC_Rabitz(
                        self.psi_coeff_t[i],
                        self.prop_chi.mol.wf.Ci(),
                        self.prop_psi.mol.MuT() + self.prop_psi.pcm.muLF,
                        self.alpha_t[i - 1])
                    self.prop_chi.propagate_one_step(self.prop_field.get_field(),
                                                     self.psi_coeff_t[i])

                    self.chi_coeff_t[i - 1] = self.prop_chi.mol.wf.Ci()

                for i in range(self.nstep):
                    if i == 0:
                        self.prop_psi.mol.wf.set_ci(self.initial_c0, 0)
                        self.psi_coeff_t[i] = self.prop_psi.mol.wf.Ci()

                    self.field_vector_t[i] = self.prop_field.propagate_field_OC_Rabitz(
                        self.prop_psi.mol.wf.Ci(),
                        self.chi_coeff_t[i],
                        self.prop_psi.mol.MuT() + self.prop_psi.pcm.muLF,
                        self.alpha_t[i])
                        # qui invece scorre
                    self.prop_psi.propagate_one_step(self.prop_field.get_field())
                    self.psi_coeff_t[i + 1] = self.prop_psi.mol.wf.Ci()  # coefficients are stored

                self.check_convergence()

                self.save.save_every_iteration(self.current_iteration, [".log", "_final_pop.dat"], np.real(self.calc_every_iteration_output()))
                self.save.save_every_n_iterations(self.current_iteration,
                                                  ["_pop_t.dat", "_field_t.dat"],
                                                  [np.real(af.population_from_wf_matrix(self.psi_coeff_t)), self.field_vector_t],
                                                  "_field_bkp.dat",
                                                  self.field_vector_t)
                self.current_iteration += 1

            else:
                self.psi_coeff_t = self.prop_psi.propagate_n_step(self.nstep,
                                                                      self.field_vector_t)

                self.save.save_every_iteration(self.current_iteration, [".log", "_final_pop.dat"], np.real(self.calc_every_iteration_output()))
                self.current_iteration += 1

        self.save.save_every_n_iterations(self.current_iteration,
                                          ["_pop_t.dat", "_field_t.dat"],
                                          [np.real(af.population_from_wf_matrix(self.psi_coeff_t)), self.field_vector_t],
                                          "_field_bkp.dat",
                                          self.field_vector_t)


    def set_initial_c0(self):
        self.initial_c0 = self.prop_psi.mol.wf.Ci()
        if self.state_to_project_on.size != self.initial_c0.size:
            sys.exit(
                "ERROR: Wrong dimension of target state: number of coefficients is different from wavefunction. ")

    def set_starting_field(self, starting_field):
        self.field_vector_t = np.copy(starting_field.get_field())


class SaveOCGenetic():


class SaveOC(Save):

    def __init__(self):
        super().__init__()

    def set_n_save(self, n_vector, n_restart):
        self.save_every_n_vector = n_vector
        self.save_n_restart = n_restart

    def save_every_iteration(self, iteration, names, vector):
        for i in range(len(names)):
            self.save_iteration_vector(iteration, names[i], vector[i])

    def save_every_n_iterations(self, iteration, names, matrix, name_restart, restart):
        for i in range(len(names)):
            if iteration % self.save_every_n_vector[i] == 0:
                self.save_iteration_matrix(names[i], iteration, matrix[i])
        if iteration % self.save_n_restart == 0:
            self.save_restart(name_restart, restart)

    def creation_save_files(self, names, headers, input_par_file):
        if input_par_file.get_oc_par()["restart"] == "false" or "norestart":
            for i in range(len(names)):
                self.create_bkp_file(names[i])
        for i in range(len(names)):
            self.print_date_header(names[i],
                                   headers[i])
        self.print_log_header(input_par_file)

    def print_log_header(self, input_par_file):
        nameT = self.folder + self.namefile + ".log"
        f = open(nameT, 'a')
        f.write("#dt: " + input_par_file.get_sys_par()["dt"] +
                " env: " + input_par_file.get_env_par()["env"] +
                " alpha: " + input_par_file.get_oc_par()["alpha"] +
                " alpha value: " + input_par_file.get_oc_par()["alpha0"] + "\n")
        if input_par_file.get_oc_par()["restart"] == "norestart":
            f.write("#WARNING: restart asked but restart file not found. Starting from guess field" + "\n")
            # input_par_file.modify_oc_restart_par("false")
        if input_par_file.get_oc_par()["restart"] != "true":
            f.write("#field parameters: \n #field: " + input_par_file.get_field_par()["field_type"] + " " +
                    "fi: " + np.array2string(input_par_file.get_field_par()["fi"]))
            if input_par_file.get_field_par()["field_type"] != "sum":
                f.write(" fi_cos: " + np.array2string(input_par_file.get_field_par()["fi_cos"]) +
                        "\n")
                f.write("# omega: " + np.array2string(input_par_file.get_field_par()["omega"]))
            elif input_par_file.get_field_par()["field_type"] != "const":
                f.write(" sigma: " + input_par_file.get_field_par()["sigma"] +
                        " omega: " + input_par_file.get_field_par()["omega"] +
                        " t0: " + input_par_file.get_field_par()["t0"])
        else:
            f.write("#Restarted from bkp field ")
        f.write("\n")
        f.write("#target state: " + input_par_file.get_oc_par()["state_to_project_on"] + "\n")
        f.write("\n")

    def print_iteration_separator(self, name):
        f = open(name, 'a')
        f.write("#new iteration \n")
        f.close()

    def save_iteration_vector(self, iteration, name, vector):
        nameT = self.folder + self.namefile + name
        f = open(nameT, 'ab')
        np.savetxt(f, np.insert(vector, 0, iteration)[None], delimiter=' ', header='', footer='', fmt='%1.8f')
        if vector.ndim > 1:
            self.print_iteration_separator(nameT)
        f.close()

    def save_iteration_matrix(self, name, iteration, matrix):
        nameT = self.folder + self.namefile + name
        f = open(nameT, 'ab')
        if matrix.ndim != 1:
            t = np.arange(matrix.shape[0])
            iteration = np.full(matrix.shape[0], iteration)
            self.print_iteration_separator(nameT)
            np.savetxt(f, np.real(np.insert(np.insert(matrix, 0, t, axis=1), 0, iteration, axis=1)),
                       delimiter=' ', header='', footer='', fmt='%1.5f')
        else:
            np.savetxt(f, np.insert(matrix, 0, iteration)[None], delimiter=' ', header='', footer='', fmt='%1.8f')
        f.close()


#    def save_J_terms(self, namefile, iteration, wavefunction, target_state, convergence, J, field_matrix, alpha_t, dt):
#        out_projection = af.projector_mean_value(wavefunction,
#                                                 target_state)/(np.dot(wavefunction, np.conj(wavefunction)))
#        out_field = af.field_J_integral(field_matrix, dt)
#        out_alpha_field = af.alpha_field_J_integral(field_matrix, alpha_t, dt)
#        f = open(namefile, 'ab')
#        out = np.real(np.hstack((iteration, convergence, J, out_projection, out_alpha_field, out_field)))
#        np.savetxt(f, out[None], delimiter=' ', header='', footer='', fmt='%1.8f')
#        f.close()


class SaveBO(Save):
    def __init__(self):
        super().__init__()

    def save_BO(self, namefile, BO_en, BO_eigvect):
        self.create_bkp_file(namefile + "_BO_energies.dat")
        self.print_date_header(namefile + "_BO_energies.dat",
                               "#Eigenvalues: \n")
        f = open(namefile + "_BO_energies.dat", 'ab')
        np.savetxt(f, BO_en, delimiter=' ', header='', footer='', fmt='%1.8f')
        f.close()
        f = open(namefile + "_BO_energies.dat", 'a')
        f.write("#Egenvectors: \n")
        f.close()
        f = open(namefile + "_BO_energies.dat", 'ab')
        np.savetxt(f, BO_eigvect, delimiter=' ', header='', footer='', fmt='%1.8f')
        f.close()


class Genetic_iterator():
    def __init__(self):

        self.genetic = GeneticParameters()
        self.prop_psi = None
        self.convergence_t = None
        self.J = None


    def initialize_propagation_attributes(self, molecule, field, env, dt, nstep, target_state):
        prop_psi = Propagator()
        prop_psi.set_propagator(dt, molecule, env)
        self.prop_psi = np.ndarray(shape=(self.genetic.n_genes), dtype=Propagator())
        for i in range(self.genetic.n_genes):
            self.prop_psi[i] = prop_psi
        self.psi_coeff_t = np.zeros([self.genetic.n_genes, nstep + 1, self.prop_psi[0].mol.wf.n_ci], dtype=complex)
        self.J = np.zeros([self.genetic.n_genes], dtype=complex)
        self.set_initial_c0(target_state)


    def set_initial_c0(self, target_state):
        self.initial_c0 = self.prop_psi[0].mol.wf.Ci()
        if target_state.size != self.initial_c0.size:
            sys.exit("ERROR: Wrong dimension of target state: number of coefficients is different from wavefunction.")



    def calc_J(self, target_state, alpha_t, dt):
        for i in range(self.genetic.n_genes):
            self.J[i] = np.real(af.projector_mean_value(self.prop_psi[i].mol.wf.Ci(), target_state) \
                 - af.alpha_field_J_integral(self.field_vector_t[i], alpha_t, dt))


    def check_convergence(self, target_state, alpha_t, dt):
        J_prev_tmp = np.copy(self.J)
        self.calc_J(target_state, alpha_t, dt)
        self.convergence_t = self.J[0] - J_prev_tmp[0]

    def get_convergence(self):
       return self.convergence_t


    def set_starting_field(self, starting_field):
        self.field_vector_t = starting_field


    def calc_field_integral(self, dt):
        integral_field = af.field_J_integral(self.field_vector_t, dt)
        return np.real(integral_field)