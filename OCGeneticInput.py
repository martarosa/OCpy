from OCInput import OCInput


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