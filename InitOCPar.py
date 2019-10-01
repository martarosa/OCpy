from read import auxiliary_functions as af
from InitPar import InitPar
from OCParameters import OCParameters
from GeneticParameters import GeneticParameters
from read.ReadNamelistGenetic import ReadNamelistGenetic


class InitOCPar(InitPar):
    def __init__(self):
        self.parameters = OCParameters()
        self.iterator_parameters = None

    def init(self, user_input):
        self.parameters.oc_iterator_name = user_input.sys.par['propagation']
        self.parameters.alpha = user_input.oc.par['alpha']
        self.parameters.alpha0 = float(user_input.oc.par['alpha0'])
        self.parameters.n_iterations = int(user_input.oc.par['n_iterations'])
        self.parameters.convergence_thr = float(user_input.oc.par['convergence_thr'])

        self.parameters.restart = user_input.oc.par['restart']

        self.parameters.nstep = int(user_input.sys.par['nstep'])
        self.parameters.dt = float(user_input.sys.par['dt'])
        self.parameters.target_state = af.normalize_vector([float(i) for i in user_input.oc.par["target_state"].split(' ')])
        self.parameters.iterator_config_file = user_input.oc.par['iterator_config_file']
        if self.parameters.oc_iterator_name == "genetic":
            self.iterator_parameters = GeneticParameters()
            genetic_input = ReadNamelistGenetic()
            if user_input.oc.par['iterator_config_file'] == 'None':
                user_input.oc.par['iterator_config_file'] = 'genetic.conf'
            genetic_input.read_file(user_input.sys.par['folder'], user_input.oc.par['iterator_config_file'])
            self.init_genetic(genetic_input)


    def init_genetic(self, genetic_input):
        self.iterator_parameters.amplitude_min = float(genetic_input.genetic.par['amplitude_min'])
        self.iterator_parameters.amplitude_max = float(genetic_input.genetic.par['amplitude_max'])
        self.iterator_parameters.n_chromosomes = int(genetic_input.genetic.par['n_chromosomes'])
        self.iterator_parameters.n_evolved_chr = int(genetic_input.genetic.par['n_evolved_chr'])