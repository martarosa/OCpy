class OCIteratorParameters():

    def __init__(self):
        self.nstep = None
        self.dt = None
        self.target_state = None
        self.alpha_t = None
        self.convergence_t =  99999
        self.J = 99999
        self.field_psi_matrix = None
        self.psi_coeff_t = None
        self.dict_out = {}
#        self.dict_restart = {}






class OCIterator():
    def __init__(self):
        self.oc_iterator_parameters = OCIteratorParameters()

    def iterate(self, current_iteration):
        pass

    def check_convergence(self):
        pass

    def calc_J(self):
        pass

    def init_output_dictionary(self):
        pass

    def init(self, oc_iterator_parameters, molecule, starting_field, pcm, alpha_t):
        pass