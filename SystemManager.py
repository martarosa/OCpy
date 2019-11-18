from field.InitFieldPar import InitFieldPar
from save.InitLogPar import InitLogPar
from molecule.InitMolecularPar import InitMolecularPar
from InitOCPar import InitOCPar
from pcm.InitPCMPar import InitPCMPar
from save.InitSavePar import InitSavePar
from read.ReadNamelistOC import ReadNamelistOC
from field.ReadFieldRestartGenetic import ReadFieldRestartGenetic
from field.ReadFieldRestartRabitz import ReadFieldRestartRabitz

from molecule.Molecule import Molecule
from field.Field import Field
from pcm.PCM import PCM, FrozenSolventPCM, DinamicPCM
from OCManager import OCManager



#The system contains Molecule, Field, PCM and OCManager objects.
# Description of the molecule, starting field and pcm are stored separately outside the OCManager. They are never modified and practically useless, there just for
# mental order and debugging

# ReadNamelistOC reads the OCinput file (no gaussian output or other data files), performs checks and store all the parameters in separate dictionaries for each section:
# sys, wf, field, env, save, oc
# Files containingmolecule and pcm informations (energies, dipoles, potential, tessere etc) are read by a separate object which knows the format (right now only ReadOutputgaussian is implemented)

# Both SystemManager, and InitPar objects (InitMolecularpar, InitiFieldpar, etc) know user_input format

# Molecule, Field, PCM, OCManager, Save objects don't know anything about user_input format.
# InitPar objects bridge the external format of user input and internal one:
# e.g. 1) InitMolecularPar.init(user_input) receives user_input as argument, read informations, open and read gaussian output files if needed
#      2) fills MoleculeParameters object which has as attributes all informations about the molecule and no methods.
#         Parameters don't know any format, are only containers of informations
#      3) Molecule.init_molecule(MoleculeParameters) initialize Molecule attributes from MoleculeParameters


class SystemManager():

    def __init__(self):

        self.mol = Molecule()
        self.starting_field = Field()
        self.pcm = PCM()
        self.oc = OCManager() # the possibility to perform a single propagation without OC is a special case of optimalControl (since this is a OC program



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
        self.mol.init_molecule(init_mol.parameters)

    def init_starting_field(self, user_input):
        init_field = InitFieldPar()
        if user_input.sys.par['propagation'] == 'genetic':
            if(user_input.field.par['field_type']) != 'genetic':
                user_input.field.par['internal_check_genetic_field'] = False
                user_input.field.par['field_type'] = 'genetic'
            init_field.read_restart = ReadFieldRestartGenetic()
        else:
            init_field.read_restart = ReadFieldRestartRabitz()
        init_field.init(user_input)
        self.starting_field.init_field(init_field.parameters)

    def init_pcm(self, user_input):
        init_pcm = InitPCMPar()
        init_pcm.init(user_input)
        if user_input.env.par['env'] == 'sol':
            self.pcm = FrozenSolventPCM()
        elif user_input.env.par['env'] == 'nanop':
            self.pcm = DinamicPCM()
        self.pcm.init_pcm(init_pcm.parameters, self.mol, self.starting_field.field[0])

    def init_optimal_control(self, user_input):
        init_oc = InitOCPar()
        init_oc.init(user_input)
        init_save = InitSavePar()
        init_save.init(user_input)
        init_log_header = InitLogPar()
        init_log_header.init(user_input)
        self.oc.init_oc(init_oc.parameters, init_oc.iterator_parameters, init_save.parameters, init_log_header.parameters, self.mol, self.starting_field, self.pcm)







