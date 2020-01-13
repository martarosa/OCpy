class OCInput:
    def __init__(self):
        self.oc_iterator_name = None
        self.alpha = None
        self.alpha0 = None
        self.n_iterations = None
        self.convergence_thr = None
        self.restart = None
        self.nstep = None
        self.dt = None
        self.target_state = None
        self.iterator_config_file = None




class OCGeneticInput(OCInput):
    def __init__(self):
        super().__init__()
        self.n_chromosomes = None
        self.n_evolved_chr = None
        self.amplitude_min = None
        self.amplitude_max = None
        self.mate = None
        self.mutate = None
        self.select = None
        self.DEAP = None


