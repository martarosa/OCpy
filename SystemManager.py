from field.SetFieldInput import SetFieldInput
from save.SetLogInput import SetLogInput
from molecule.SetMoleculeInput import SetMoleculeInput
from SetOCInput import SetOCInput
from pcm.SetPCMInput import SetPCMInput
from save.SetSaveInput import SetSaveInput
from read.ReadNamelistOC import ReadNamelistOC
from field.ReadFieldRestartGenetic import ReadFieldRestartGenetic
from field.ReadFieldRestartRabitz import ReadFieldRestartRabitz

from molecule.Molecule import Molecule
from field.Field import Field
from pcm.ABCPCM import ABCPCM
from pcm.DinamicPCM import DinamicPCM
from pcm.FrozenSolventPCM import FrozenSolventPCM
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
        self.pcm = None #ABCPCM()
        self.oc = OCManager() # the possibility to perform a single propagation without OC is a special case of optimalControl (since this is a OC program



    def init_system(self, folder, name_file):

        user_input = ReadNamelistOC() #tmp, after reading everithing apart system manager vanishes
        user_input.read_file(folder, name_file)

        self.init_molecule(user_input)
        self.init_starting_field(user_input)
        if user_input.env.section_dictionary['env'] != "vac":
            self.init_pcm(user_input)
        self.init_optimal_control(user_input)


    def init_molecule(self, user_input):
        init_mol = SetMoleculeInput()
        init_mol.set(user_input)
        self.mol.init_molecule(init_mol.input_parameters)

    def init_starting_field(self, user_input):
        init_field = SetFieldInput()
        if user_input.sys.section_dictionary['oc_algorithm'] == 'genetic':
            init_field.read_restart = ReadFieldRestartGenetic()
        else:
            init_field.read_restart = ReadFieldRestartRabitz()
        init_field.set(user_input)
        self.starting_field.init_field(init_field.input_parameters)

    def init_pcm(self, user_input):
        init_pcm = SetPCMInput()
        init_pcm.set(user_input)
        if user_input.env.section_dictionary['env'] == 'sol':
            self.pcm = FrozenSolventPCM()
        elif user_input.env.section_dictionary['env'] == 'nanop':
            self.pcm = DinamicPCM()
        self.pcm.init_pcm(init_pcm.input_parameters, self.mol, self.starting_field.field[0])

    def init_optimal_control(self, user_input):
        init_oc = SetOCInput()
        init_oc.set(user_input)
        init_save = SetSaveInput()
        init_save.set(user_input)
        init_log_header = SetLogInput()
        init_log_header.set(user_input)
        self.oc.init_oc(init_oc.input_parameters, init_save.input_parameters, init_log_header.input_parameters, self.mol, self.starting_field, self.pcm)







