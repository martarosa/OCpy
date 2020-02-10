from abc import ABCMeta, abstractmethod



class OCIteratorParameters():
    def __init__(self):
        self.target_state = None
        self.alpha_t = None


class SimulationParameters():
    def __init__(self):
        self.nstep = None
        self.dt = None


class ABCOCIterator(metaclass=ABCMeta):
    def __init__(self):
        self.par = OCIteratorParameters()
        self.simulation_par = SimulationParameters()

        self.convergence_t = None
        self.J = None
        self.field_psi_matrix = None
        self.psi_coeff_t = None
        self.dict_out = {}


    @abstractmethod
    def iterate(self, current_iteration):
        pass

    @abstractmethod
    def check_convergence(self):
        pass

    @abstractmethod
    def calc_J(self):
        pass

    @abstractmethod
    def init(self, molecule, starting_field, pcm, alpha_t, oc_input, iterator_config_input):
        pass

    @abstractmethod
    def init_output_dictionary(self):
        pass


    @abstractmethod
    def get_restart(self):
        pass
