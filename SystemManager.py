from read_and_set.read.ReadNamelistOC import ReadNamelistOC
from read_and_set.read.ReadNamelistGenetic import ReadNamelistGenetic

from read_and_set.set.SetFieldInput import SetFieldInput
from read_and_set.read.ReadFieldRestartGenetic import ReadFieldRestartGenetic
from read_and_set.read.ReadFieldRestartRabitz import ReadFieldRestartRabitz
from read_and_set.set.SetMoleculeInput import SetMoleculeInput
from read_and_set.set.SetMediumInput import SetMediumInput
from read_and_set.set.SetOCInput import SetOCInput
from read_and_set.set.SetGeneticOCInput import SetGeneticOCInput
from read_and_set.set.SetNoConfOCInput import SetNoConfOCInput
from read_and_set.set.SetSaveInput import SetSaveInput
from read_and_set.set.SetLogInput import SetLogInput

from molecule.Molecule import Molecule
from field.Field import Field
from OCManager import OCManager
from dictionaries import Dictionaries as dict

#The system contains Molecule, Field, PCM and OCManager objects.
# Description of the molecule, starting field and medium are stored separately outside the OCManager. They are never modified and practically useless, there just for
# mental order and debugging

# ReadNamelistOC reads the OCinput file (no gaussian output or other data files), performs checks and store all the parameters in separate dictionaries for each section:
# sys, wf, field, medium, save, oc
# Files containingmolecule and medium informations (energies, dipoles, potential, tessere etc) are read by a separate object which knows the format (right now only ReadOutputgaussian is implemented)

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
        self.medium = None #ABCPCM()
        self.oc = OCManager() # the possibility to perform a single propagation without OC is a special case of optimalControl (since this is a OC program




    def init_system(self, folder, name_file):
        user_input = ReadNamelistOC() #tmp, after reading everithing apart system manager vanishes
        user_input.read_file(folder, name_file)
        self.init_molecule(user_input)
        self.init_starting_field(user_input)
        self.init_medium(user_input)
        self.init_optimal_control(user_input)


    def init_molecule(self, user_input):
        init_mol = SetMoleculeInput()
        init_mol.set(user_input)
        self.mol.init_molecule(init_mol.input_parameters)

    def init_starting_field(self, user_input):
        set_field = SetFieldInput()
        if user_input.sys.section_dictionary['oc_algorithm'] == 'genetic':
            set_field.read_restart = ReadFieldRestartGenetic()
        else:
            set_field.read_restart = ReadFieldRestartRabitz()
        set_field.set(user_input)
        self.starting_field.init_field(set_field.input_parameters)


    def init_medium(self, user_input):
        print(user_input.medium.section_dictionary['medium'])
        set_medium = SetMediumInput()
        self.medium = dict.MediumDict[user_input.medium.section_dictionary['medium']]()
        set_medium.set(user_input)
        self.medium.init_medium(set_medium.input_parameters, self.mol, self.starting_field.field)


    def init_optimal_control(self, user_input):
        set_oc = SetOCInput()
        set_oc.set(user_input)
        if user_input.sys.section_dictionary['oc_algorithm'] == 'genetic':
            iterator_config_input = ReadNamelistGenetic()
            iterator_config_input.read_file(user_input.sys.section_dictionary['folder'],
                                            user_input.oc.section_dictionary['iterator_config_file'])
            set_iterator_config = SetGeneticOCInput()
            set_iterator_config.set(iterator_config_input)
        else:
            iterator_config_input =  None
            set_iterator_config = SetNoConfOCInput()

        set_save = SetSaveInput()
        set_save.set(user_input)
        set_log_header = SetLogInput()
        set_log_header.set(user_input)
        if user_input.sys.section_dictionary['oc_algorithm'] == 'genetic':
            set_log_header.set_conf_log(iterator_config_input.string_file_config)


        self.oc.init_oc(set_oc.input_parameters,
                        set_iterator_config.input_parameters,
                        set_save.input_parameters,
                        set_log_header.input_parameters,
                        self.mol,
                        self.starting_field,
                        self.medium)







