import numpy as np
import time
import os

class SaveTools():

    def create_bkp_file(self, name_file):
        if os.path.isfile(name_file):
            n_bkp = 1
            bkp_name = name_file + "." + str(n_bkp)
            while os.path.isfile(bkp_name):
                n_bkp += 1
                bkp_name = name_file + "." + str(n_bkp)
            for i in range(n_bkp, 1, -1):
                os.rename(name_file + "." + str(i - 1), name_file + "." + str(i))
            os.rename(name_file, name_file + ".1")

    def print_date_header(self, name_file, header):
        try:
            f = open(name_file, 'a')
        except IOError:
            f = open(name_file, 'w+')
        f.write("#Run started at " + time.strftime("%H:%M:%S") +" " + time.strftime("%d/%m/%Y") +"\n")
        f.write(header)
        f.close()

    def print_iteration_separator(self, name_file):
        f = open(name_file, 'a')
        f.write("#new iteration \n")
        f.close()

    def save_iteration_matrix(self, name_file, iteration, matrix):
        f = open(name_file, 'ab')
        if matrix.ndim > 1:
            t = np.arange(matrix.shape[0])
            iteration = np.full(matrix.shape[0], iteration)
            self.print_iteration_separator(name_file)
            np.savetxt(f, np.real(np.insert(np.insert(matrix, 0, t, axis=1), 0, iteration, axis=1)),
                       delimiter=' ', header='', footer='', fmt='%1.5f')
        else:
            np.savetxt(f, np.insert(matrix, 0, iteration)[None], delimiter=' ', header='', footer='', fmt='%1.8f')
        f.close()

    def print_log_header(self, name_file, log_header_par):
        #check log extension
        f = open(name_file+".log", 'a')
        f.write("#dt: "+ log_header_par.dt +
                " env: " + log_header_par.env +
                " alpha: " + log_header_par.alpha +
                " alpha value: " + log_header_par.alpha0 + "\n")
        if log_header_par.restart == "norestart_found":
            f.write("#WARNING: restart asked but restart file not found. Starting from guess field" + "\n")
            #input_par_file.modify_oc_restart_par("false")
        if log_header_par.restart == "norestart_found":
            f.write("#field parameters: \n #field: " + log_header_par.field_type + " fi: " + log_header_par.fi)
            if log_header_par.field_type == "sum":
                f.write(" fi_cos: " + log_header_par.fi_cos + "\n")
                f.write("# omega: " + log_header_par.omega)
            elif log_header_par.field_type != "const":
                f.write(" sigma: " + log_header_par.sigma +
                        " omega: " + log_header_par.omega +
                        " t0: " + log_header_par.t0)
        else:
            f.write("#Restarted from bkp field: " + name_file)
        f.write("\n")
        f.write("#target state: " + log_header_par.target_state +"\n")
        f.write("\n")

    def creation_save_files(self, folder_name, save_file, restart):
        name = folder_name + save_file.name
        if restart == "false" or restart == "norestart_found" or restart == 'restart_from_different_name_field':
            self.create_bkp_file(name)
        self.print_date_header(name, save_file.header)


    def save_every_n_iterations(self, iteration, folder_name, save_file):
        name = folder_name + save_file.name
        if iteration % save_file.nstep == 0:
            self.save_iteration_matrix(name, iteration, save_file.out())

    def save_restart(self, iteration, folder_name, save_restart):
        name = folder_name + save_restart.name
        if iteration % save_restart.nstep == 0:
            np.savetxt(name, save_restart.out(), delimiter =' ', header ='', footer ='')