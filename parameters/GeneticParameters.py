class GeneticParameters():
    def __init__(self):
        self.genetic_algorithm = None
        self.n_chromosomes = None
        self.n_selected_chr = None

        self.amplitude_lim = None
        self.n_amplitudes = None
        self.omegas_matrix = None
        self.mate = None
        self.mate_probability = None

        self.mutate = None
        self.mutate_probability = None
        self.mutate_mu = None
        self.mutate_starting_sigma = None
        self.select = None