class OCIterator():
    def __init__(self):
        self.nstep = None
        self.dt = None

        self.field_psi_matrix = None
        self.psi_coeff_t = None

        self.alpha_t = None
        self.convergence_t =  99999
        self.J = 99999
        self.dict_out = {}
        self.dict_restart = {}


    def iterate(self, current_iteration):
        pass

    def check_convergence(self):
        pass

    def calc_J(self):
        pass

    def init_output_dictionary(self):
        pass

    def init_oc_iterator(self, molecule, starting_field, env, oc_parameters, alpha_t):
        pass