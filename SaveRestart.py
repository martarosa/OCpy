class SaveRestart():
    def __init__(self, name, nstep, oc_iterator):
        self.name = name
        self.nstep = nstep
        self.out = oc_iterator.get_restart