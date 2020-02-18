class SaveRestart():
    def __init__(self, name, restart_step, oc_iterator):
        self.name = name
        self.restart_step = restart_step
        self.out = oc_iterator.get_restart