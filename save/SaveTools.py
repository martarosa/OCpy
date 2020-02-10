import time
import os

import numpy as np



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


    def creation_save_files(self, folder_name, save_file, append):
        name = folder_name + save_file.name
        if append == "true":
            if not os.path.isfile(name):
                print("WARNING: file " + save_file.name + " not found in folder. Writing in new file")
        else:
            self.create_bkp_file(name)
        self.print_date_header(name, save_file.header)

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

    def save_every_n_iterations(self, iteration, folder_name, save_file):
        name = folder_name + save_file.name
        if iteration % save_file.nstep == 0:
            self.save_iteration_matrix(name, iteration, save_file.out())

    def save_restart(self, iteration, folder_name, save_restart):
        name = folder_name + save_restart.name
        if iteration % save_restart.nstep == 0:
            np.savetxt(name, save_restart.out(), delimiter =' ', header ='', footer ='')

    def save_3D_matrix(self, Mijn, name_file):
        f = open(name_file, 'w+')
        for i in range(Mijn.shape[0]):
            for j in range(Mijn.shape[1]):
                f.write(str(i)+" "+str(j)+"\n")
                np.savetxt(f, Mijn[i,j])
        f.close()
