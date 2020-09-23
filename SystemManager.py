import dictionaries.MediumDictionaries
import dictionaries.OCDictionaries
import dictionaries.FieldDictionary

from OCManager import OCManager
from field.Field import Field
from molecule.Molecule import Molecule


from read_and_set.read.input_sections.ReadInput import ReadInput
from read_and_set.set.SetFieldInput import SetFieldInput
from read_and_set.set.SetLogInput import SetLogInput
from read_and_set.set.SetMediumInput import SetMediumInput
from read_and_set.set.SetMoleculeInput import SetMoleculeInput
from read_and_set.set.SetOCInput import SetOCInput
from read_and_set.set.SetSaveInput import SetSaveInput


#The system contains Molecule, Field, Medium and OCManager objects.
# Molecule, Field and Medium are both here and inside OCManager. These ones are never modified and here only for debugging


# In init_system method ReadInput() reads the input file performs checks and store all the parameters in separate
#  dictionaries for each of the input sections: sys, wf, field, medium, save, oc

# Then Set*() methods extract all the information needed for the calulation and fill correspondinf data structures,
# used to finally initialize internal objects. Read and Set methods know the input keys, flags and structure, after them
# the information is sitributed stored in internal classes. More informations on read and set are in read_adn_set README


class SystemManager():

    def __init__(self):

        self.mol = Molecule()
        self.starting_field = Field()
        self.medium = None
        self.oc = OCManager()



    def init_system(self, folder, name_file):
        user_input = ReadInput()
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
        set_field.read_restart = dictionaries.FieldDictionary.FieldRestartDict[user_input.sys.section_dictionary['oc_algorithm']]()
        set_field.set(user_input)
        self.starting_field.init_field(set_field.input_parameters)

    def init_medium(self, user_input):
        print(user_input.medium.section_dictionary['medium'])
        set_medium = SetMediumInput()
        self.medium = dictionaries.MediumDictionaries.MediumDict[user_input.medium.section_dictionary['medium']]()
        set_medium.set(user_input)
        self.medium.init_medium(set_medium.input_parameters, self.mol, self.starting_field.field)



    def init_optimal_control(self, user_input):
        set_oc = SetOCInput()
        set_oc.set(user_input)
        OCConf = dictionaries.OCDictionaries.OCAlgorithmConfig[user_input.sys.section_dictionary['oc_algorithm']]()
        set_OCConf = dictionaries.OCDictionaries.OCAlgorithmSet[user_input.sys.section_dictionary['oc_algorithm']]()
        if not user_input.sys.section_dictionary['oc_algorithm'] == "none":
            OCConf.read_file(user_input.sys.section_dictionary['folder'],
                                            user_input.oc.section_dictionary['conf_file'])
        set_OCConf.set(OCConf)
        set_save = SetSaveInput()
        set_save.set(user_input)
        set_log_header = SetLogInput()
        set_log_header.set(user_input)
        if user_input.sys.section_dictionary['oc_algorithm'] != 'none':
            set_log_header.set_conf_log(OCConf.conf_str)


        self.oc.init_oc(set_oc.input_parameters,
                        set_OCConf.input_parameters,
                        set_save.input_parameters,
                        set_log_header.input_parameters,
                        self.mol,
                        self.starting_field,
                        self.medium)







