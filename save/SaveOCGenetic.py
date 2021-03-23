import numpy as np

from save.SaveFile import SaveFile
from save.ABCSave import ABCSave
from parameters.SaveParameters import SaveParameters
from save.SaveRestart import SaveRestart
from save.SaveTools import SaveTools

class SaveOCGenetic(ABCSave):
    def __init__(self):
        super().__init__()
        self.par = SaveParameters()
        self.save_files = []
        self.restart_file = None

        self.save_tools = SaveTools()


    def init_save_file_list(self, oc_iterator, restart_step):
        log = SaveFile(".log",
                       "#fields: n_iteration, convergence, J, projection on target state, alpha*integral(field^2)\n",
                       1,
                       'log_file',
                       oc_iterator)

        final_pop = SaveFile("_final_pop.dat",
                             "#Final states populations \n #fields: n_iteration, final states population \n",
                             1,
                            'final_pop',
                            oc_iterator)

        pop_t = SaveFile("_pop_t.dat",
                         "#pop(t) \n #fields: n_iteration, nstep, states population(t)  \n",
                         restart_step,
                         'pop_t',
                         oc_iterator)

        field_t = SaveFile("_field_t.dat",
                           "#field(t) \n#fields: n_iteration, nstep, time, field(t) x, y, z  \n",
                           restart_step,
                           'field_t',
                           oc_iterator)

        field_ampl = SaveFile("_field_ampl.dat",
                              "#field_amplitudes \n #field shape = a0 + sum(ai*sin(wi t)) \n "
                              +"#omegas = "+np.array_repr(oc_iterator.genetic_par.omegas_matrix[:,0]).replace('\n', ''),
                              restart_step,
                              'field_ampl',
                              oc_iterator)

        self.save_files = [log, final_pop, pop_t, field_t, field_ampl]
        self.restart_file = SaveRestart("_field_bkp.dat", restart_step, oc_iterator)



    def init_save(self, save_input, log_input, oc_iterator):
        self.par.folder = save_input.folder
        self.par.filename = save_input.name
        self.par.restart_calculation = save_input.restart_calculation
        self.par.append = save_input.append
        self.init_save_file_list(oc_iterator, save_input.restart_step)


        for i in np.arange(len(self.save_files)):
            print("save: " + str(i))
            self.save_tools.creation_save_files(self.par.folder + self.par.filename, self.save_files[i], self.par.append)
        self.print_log_header(log_input)



    def save(self, iteration):
        for i in np.arange(len(self.save_files)):
            self.save_tools.save_every_n_iterations(iteration, self.par.folder + self.par.filename, self.save_files[i])
        self.save_tools.save_restart(iteration, self.par.folder + self.par.filename, self.restart_file)
