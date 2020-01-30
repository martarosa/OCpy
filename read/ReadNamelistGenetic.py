import configparser
import sys
from read import NameListSections as sec
from read import ReadNamelist as readnamelist



class ReadNamelistGenetic(readnamelist.ReadNamelist):
    def __init__(self):
        super().__init__()
        self.n_sections = None
        self.default_dict = None
        self.sections = ["GENETIC", "MATE", "MUTATE", "SELECT"]
        self.genetic = sec.SectionGenetic()
        self.mate = sec.SectionMate()
        self.mutate = sec.SectionMutate()
        self.select = sec.SectionSelect()


    def set_default_dict(self):
        default_dict = ({'GENETIC': {
            'chromosomes': '120',
            'n_evolver_chr': '20',
            'genetic_algorithm' : 'sequential',
            'amplitude_lim': '0.05'},
            'MATE': {
                'mate' : 'cxUniform',
                'mate_probability': '1'
            },
                'MUTATE': {
                'mutate':'mutGaussian',
                'mutate_probability': '0.2',
                'n_mutate': '2',
                'starting_sigma': '0.01',
                'eta_thr': '0.6',
                'q': '0.9'
            },
            'SELECT':{
                'select':'selBest'
            }
        })



    def read_file(self, folder, namefile):
        user_input = configparser.ConfigParser()
        #first read without default values for checks on consistency
        user_input.read(folder + namefile)
        self.check_input_sections_names(user_input)
        self.check_input_consistency(user_input)
        #after checks on consistency we read again with default
        user_input.read_dict(self.default_dict)
        self.n_sections = len(user_input.sections())
        user_input.read(folder + namefile)
        self.set_nml_sections(user_input)


    def set_nml_sections(self, user_input):
        self.genetic.section_dictionary = dict(user_input.items("GENETIC"))
        self.genetic.lowering_case()
        self.genetic.check_allowed_values()
        self.mate.section_dictionary = dict(user_input.items("MATE"))
        self.mate.lowering_case()
        self.mate.check_allowed_values()
        self.mutate.section_dictionary = dict(user_input.items("MUTATE"))
        self.mutate.lowering_case()
        self.mutate.check_allowed_values()
        self.select.section_dictionary = dict(user_input.items("SELETC"))
        self.select.lowering_case()
        self.select.check_allowed_values()



    def check_input_consistency(self, user_input):
        self.check_genetic_nml_consistency(user_input)
        self.check_mate_nml_consistency(user_input)
        self.check_mutate_nml_consistency(user_input)
        self.check_select_nml_consistency(user_input)



    def check_genetic_nml_consistency(self, user_input):
        pass

    def check_mate_nml_consistency(self, user_input):
        pass

    def check_mutate_nml_consistency(self, user_input):
        #if propagation OPTIMALC namelist mus be empty
        if (user_input['GENETIC']['genetic_algorithm'].lower() != 'mixed'):

            if len(user_input['OPTIMALC'].keys()) != 0:
                sys.exit("Error. Propagation without optimal control. OPTIMALC namelist should be empty")

    def check_select_nml_consistency(self, user_input):
        pass












class ReadOCConfigurationNamelistGenetic():
    def __init__(self):
        self.n_sections = None
        self.genetic = sec.SectionGenetic()
        self.mate = sec.SectionMate()
        self.mutate = sec.SectionMutate()
        self.select = sec.SectionSelect()

    def read_file(self, folder, namefile):
        user_input = configparser.ConfigParser()
        user_input.read_dict({'GENETIC': {
                                  'amplitude_min': '0.01',
                                  'amplitude_max': '0.01',
                                  'n_chromosomes': '1',
                                  'n_evolved_chr': 1,
                                  'mate': 'DEAP_cxUniform',
                                  'mutate': 'DEAP_mutGaussian',
                                  'select':'DEAP_selBest',
                                  'deap': 'true'}
                              })
        user_input.read(folder + namefile)
        self.check_input_sections(user_input)
        self.set_dictionaries(user_input)




    def read_file_new(self, folder, namefile):
        user_input = configparser.ConfigParser()
        default_dict = ({'GENETIC': {
            'chromosomes': '120',
            'n_evolver_chr': '20',
            'genetic_algorithm' : 'sequential',
            'amplitude_lim': '0.05'},
            'MATE': {
                'mate' : 'cxUniform',
                'mate_probability': '1'
            },
                'MUTATE': {
                'mutate':'mutGaussian',
                'mutate_probability': '0.2',
                'n_mutate': '2',
                'starting_sigma': '0.01',
                'eta_thr': '0.6',
                'q': '0.9'
            },
            'SELECT':{
                'select':'selBest'
            }
        })

        #first read without default values for checks on consistency
        user_input.read(folder + namefile)
        self.check_input_sections(user_input)
        self.check_input_consistency(user_input)
        #after checks on consistency we read again with default
        user_input.read_dict(default_dict)
        self.n_sections = len(user_input.sections())
        user_input.read(folder + namefile)
        self.set_dictionaries(user_input, default_dict)



    def set_dictionaries(self, user_input):
        self.genetic.section_dictionary = dict(user_input.items("GENETIC"))
        self.genetic.lowering_case()
        self.genetic.check_allowed_values()


    def check_input_sections(self, user_input):
        if len(user_input.sections()) != 1:
            sys.exit("Error. Sections names are wrong in input file")



    def check_input_consistency(self, user_input):
        self.check_genetic_nml_consistency(user_input)
        self.check_mate_nml_consistency(user_input)
        self.check_mutate_nml_consistency(user_input)
        self.check_select_nml_consistency(user_input)


    def check_genetic_nml_consistency(self, user_input):
        pass

    def check_mate_nml_consistency(self, user_input):
        pass

    def check_mutate_nml_consistency(self, user_input):
        #if propagation OPTIMALC namelist mus be empty
        if (user_input['GENETIC']['genetic_algorithm'].lower() != 'mixed'):

            if len(user_input['OPTIMALC'].keys()) != 0:
                sys.exit("Error. Propagation without optimal control. OPTIMALC namelist should be empty")

    def check_select_nml_consistency(self, user_input):
        pass



