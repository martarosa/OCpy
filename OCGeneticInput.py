from OCInput import OCInput


class OCGeneticInput(OCInput):
    def __init__(self):
        self.n_chromosomes = None
        self.n_selected_chr = None

        self.amplitude_lim = None
        self.mate = None
        self.mate_probability = None

        self.mutate = None
        self.mutate_probability = None
        self.n_mutate = None
        self.mutate_starting_sigma = None
        self.eta_thr = None
        self.q = None
        self.select = None