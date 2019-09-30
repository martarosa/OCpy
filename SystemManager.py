import numpy as np
import sys

import os.path



from ReadOutputGaussian import ReadOutputGaussian
from ReadNamelistOC import ReadNamelistOC
from ReadFieldRestart import ReadFieldRestartGenetic, ReadFieldRestartRabitz, ReadFieldRestart

from Molecule import Molecule
from Field import Field
from PCM import PCM, FrozenSolventPCM, DinamicPCM
from OCManager import OCManager

from MoleculeParameters import MoleculeParameters
from PCMParameters import PCMParameters
from LogHeaderParameters import LogHeaderParameters
from SaveParameters import SaveParameters
from FieldParameters import FieldParameters
from OCParameters import OCParameters

import auxiliary_functions as af
from OCParameters import GeneticParameters, InitGeneticPar


class InitPar():
    def __init__(self):
        self.parameters = MoleculeParameters()


class InitMolecularPar():
    def __init__(self):
        self.mol_parameters = MoleculeParameters()


    def init(self, user_input):
        read_output = ReadOutputGaussian()

        self.mol_parameters.wf_ci = read_output.read_ci0(user_input.sys.par['folder'] +
                                                         user_input.wf.par['name_ci'])

        self.mol_parameters.en_ci = read_output.read_en_ci0(user_input.sys.par['folder'] +
                                                            user_input.wf.par['name_ei'])
        self.mol_parameters.muT = read_output.read_muT(user_input.sys.par['folder'] +
                                                       user_input.wf.par['name_mut'], self.mol_parameters.wf_ci.size)
        if user_input.env.par['env'] != 'vac':
            self.mol_parameters.Vijn = read_output.read_V(user_input.sys.par['folder'] +
                                                          user_input.env.par['name_vij'],
                                                          self.mol_parameters.wf_ci.size)

class InitFieldPar():
    def __init__(self):
        self.field_parameters = FieldParameters()
        self.read_restart = ReadFieldRestart()

    def init(self, user_input):
        if (user_input.oc.par['restart'] == 'false'):
            self.init_no_restart(user_input)
        else:
            self.init_restart(user_input)



    def init_no_restart(self, user_input):
        read_output = ReadOutputGaussian()
        self.field_parameters.dt = float(user_input.sys.par['dt'])
        self.field_parameters.nstep = int(user_input.sys.par['nstep'])
        self.field_parameters.field_type = user_input.field.par['field_type']
        self.field_parameters.fi = user_input.field.par['fi']
        self.field_parameters.fi_cos = user_input.field.par['fi_cos']
        self.field_parameters.sigma = float(user_input.field.par['sigma'])
        self.field_parameters.omega = user_input.field.par['omega']
        self.field_parameters.t0 = float(user_input.field.par['t0'])
        self.field_parameters.omega_max = read_output.read_en_ci0(user_input.sys.par['folder'] +
                                                            user_input.wf.par['name_ei'])[-1]



    def init_restart(self, user_input):
        if os.path.isfile(user_input.sys.par['folder'] + user_input.field.par['name_field_file']):
            if(user_input.sys.par['folder'] + user_input.field.par['name_field_file'] != user_input.sys.par['folder'] + user_input.sys.par['name'] + '_field_bkp.dat'):
                user_input.oc.par['restart'] = 'restart_from_different_name_field'
            self.read_restart.read_file(user_input.sys.par['folder'], user_input.field.par['name_field_file'])
            self.field_parameters = self.read_restart.field_par
        elif os.path.isfile(user_input.sys.par['folder'] + user_input.sys.par['name'] + '_field_bkp.dat'):
            self.read_restart.read_file(user_input.sys.par['folder'], user_input.sys.par['name'] + '_field_bkp.dat')
            self.field_parameters = self.read_restart.field_par
        else:
            user_input.oc.par['restart'] = 'norestart_found'
            self.init_no_restart(user_input)




class InitPCMPar():
    def __init__(self):
        self.pcm_parameters = PCMParameters()

    def init(self, user_input):
        read_output = ReadOutputGaussian()
        self.pcm_parameters.cavity = read_output.read_cavity_tesserae(user_input.sys.par['folder'] + user_input.env.par['name_file_cavity'])
        self.pcm_parameters.env = user_input.env.par['env']
        if self.pcm_parameters.env == 'sol':
             self.pcm_parameters.Qnn_reactionfield = read_output.read_Q_matrix(user_input.sys.par['folder'] +
                                                                               user_input.env.par['name_q_tdplas'])
             self.pcm_parameters.Qnn_localfield = read_output.read_Q_matrix(user_input.sys.par['folder'] +
                                                                            user_input.env.par['name_q_local_field'])


class InitSavePar():
    def __init__(self):
        self.save_parameters = SaveParameters()

    def init(self, user_input):
        self.save_parameters.folder = user_input.sys.par['folder']
        self.save_parameters.name = user_input.sys.par['name']
        self.save_parameters.restart_calculation = user_input.oc.par['restart']
        self.save_parameters.restart_step = int(user_input.save.par['restart_step'])


class InitLogPar():
    def __init__(self):
        self.log_header_parameters = LogHeaderParameters()

    def init(self, user_input):
        self.log_header_parameters.dt = user_input.sys.par['dt']
        self.log_header_parameters.env = user_input.env.par['env']

        self.log_header_parameters.restart = user_input.oc.par['restart']
        self.log_header_parameters.target_state = user_input.oc.par["target_state"]
        self.log_header_parameters.alpha = user_input.oc.par['alpha']
        self.log_header_parameters.alpha0 = user_input.oc.par['alpha0']

        self.log_header_parameters.field_type = user_input.field.par['field_type']
        self.log_header_parameters.fi = np.array2string(user_input.field.par['fi'])
        self.log_header_parameters.fi_cos = np.array2string(user_input.field.par['fi_cos'])
        self.log_header_parameters.sigma = user_input.field.par['sigma']
        self.log_header_parameters.omega = np.array2string(user_input.field.par['omega'])
        self.log_header_parameters.t0 = user_input.field.par['t0']


class InitOCPar():
    def __init__(self):
        self.oc_parameters = OCParameters()

    def init(self, user_input):
        self.oc_parameters.oc_iterator_name = user_input.sys.par['propagation']
        self.oc_parameters.alpha = user_input.oc.par['alpha']
        self.oc_parameters.alpha0 = float(user_input.oc.par['alpha0'])
        self.oc_parameters.n_iterations = int(user_input.oc.par['n_iterations'])
        self.oc_parameters.convergence_thr = float(user_input.oc.par['convergence_thr'])

        self.oc_parameters.restart = user_input.oc.par['restart']

        self.oc_parameters.nstep = int(user_input.sys.par['nstep'])
        self.oc_parameters.dt = float(user_input.sys.par['dt'])
        self.oc_parameters.target_state = af.normalize_vector([float(i) for i in user_input.oc.par["target_state"].split(' ')])
        print(self.oc_parameters.target_state)
        if self.oc_parameters.oc_iterator_name == "genetic":
            self.oc_parameters.iterator_parameters = GeneticParameters()
            init_genetic = InitGeneticPar()
            init_genetic.init(user_input)



class SystemManager():

    def __init__(self):
        self.mol = Molecule()
        self.starting_field = Field()
        self.pcm = PCM()
        self.oc = OCManager() #l'alternativa di sola propagazione sta dentro a OC perchè questo è un programma per l'OC

    def init_system(self, folder, name_file):
        user_input = ReadNamelistOC()
        user_input.read_file(folder, name_file)
        self.init_molecule(user_input)
        self.init_starting_field(user_input)
        if user_input.env.par['env'] != "vac":
            self.init_pcm(user_input)
        self.init_optimal_control(user_input)


    def init_molecule(self, user_input):
        init_mol = InitMolecularPar()
        init_mol.init(user_input)
        self.mol.init_molecule(init_mol.mol_parameters)

    def init_starting_field(self, user_input):
        init_field = InitFieldPar()
        if user_input.sys.par['propagation'] == 'genetic':
            init_field.read_restart = ReadFieldRestartGenetic()
        else:
            init_field.read_restart = ReadFieldRestartRabitz()
        init_field.init(user_input)
        self.starting_field.init_field(init_field.field_parameters)


    def init_pcm(self, user_input):
        init_pcm = InitPCMPar()
        init_pcm.init(user_input)
        if user_input.env.par['env'] == 'sol':
            self.pcm = FrozenSolventPCM()
        elif user_input.env.par['env'] == 'nanop':
            self.pcm = DinamicPCM()
        self.pcm.init_pcm(init_pcm.pcm_parameters, self.mol, self.starting_field.field[0])


    def init_optimal_control(self, user_input):
        init_oc = InitOCPar()
        init_oc.init(user_input)
        init_save = InitSavePar()
        init_save.init(user_input)
        init_log_header = InitLogPar()
        init_log_header.init(user_input)
        self.oc.init_oc(init_oc.oc_parameters, init_save.save_parameters, init_log_header.log_header_parameters, self.mol, self.starting_field, self.pcm)







