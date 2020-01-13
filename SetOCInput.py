from read import auxiliary_functions as af
from SetInput import SetInput
from OCInput import OCInput, OCGeneticInput

from read.ReadOCConfigurationNamelistGenetic import ReadOCConfigurationNamelistGenetic


class SetOCInput(SetInput):
    def __init__(self):
        self.input_parameters = OCInput()
        #self.OC_iterator_input = None


    def set(self, user_input):
        if self.input_parameters.oc_iterator_name == "genetic":
            self.OC_iterator_input = OCGeneticInput()
        else:
            self.input_parameters = OCInput()

        self.set_common_oc_parameters(user_input)

        self.read_config_file(self.input_parameters.oc_iterator_name,
                              user_input.oc.section_dictionary['iterator_config_file'],
                                user_input.sys.section_dictionary['folder'])




            #self.init_genetic(genetic_input)


    def set_common_oc_parameters(self, user_input):
        self.input_parameters.oc_iterator_name = user_input.sys.section_dictionary['propagation']
        self.input_parameters.alpha = user_input.oc.section_dictionary['alpha']
        self.input_parameters.alpha0 = float(user_input.oc.section_dictionary['alpha0'])
        self.input_parameters.n_iterations = int(user_input.oc.section_dictionary['n_iterations'])
        self.input_parameters.convergence_thr = float(user_input.oc.section_dictionary['convergence_thr'])

        self.input_parameters.restart = user_input.oc.section_dictionary['restart']

        self.input_parameters.nstep = int(user_input.sys.section_dictionary['nstep'])
        self.input_parameters.dt = float(user_input.sys.section_dictionary['dt'])
        self.input_parameters.target_state = af.normalize_vector([float(i) for i in user_input.oc.section_dictionary["target_state"].split(' ')])
        self.input_parameters.iterator_config_file = user_input.oc.section_dictionary['iterator_config_file']



    def read_config_file(self, oc_iterator_name, config_name, folder):
        if oc_iterator_name == "genetic":
            genetic_input = ReadOCConfigurationNamelistGenetic()
            if config_name == 'None':
                config_name = 'genetic.conf'
            genetic_input.read_file(folder, config_name)
            self.init_genetic(genetic_input)
        else:
            pass #rabitz ed eulero non fanno niente

    def init_genetic(self, genetic_input):
        self.OC_iterator_input.amplitude_min = float(genetic_input.genetic.section_dictionary['amplitude_min'])
        self.OC_iterator_input.amplitude_max = float(genetic_input.genetic.section_dictionary['amplitude_max'])
        self.OC_iterator_input.n_chromosomes = int(genetic_input.genetic.section_dictionary['n_chromosomes'])
        self.OC_iterator_input.n_evolved_chr = int(genetic_input.genetic.section_dictionary['n_evolved_chr'])
        self.OC_iterator_input.DEAP = genetic_input.genetic.section_dictionary['deap']
        self.OC_iterator_input.mate = genetic_input.genetic.section_dictionary['mate']
        self.OC_iterator_input.mutate = genetic_input.genetic.section_dictionary['mutate']
        self.OC_iterator_input.select = genetic_input.genetic.section_dictionary['select']
