from propagator.PropagatorTerms import PropagatorTerms

# in propagator_terms there are all possible propagation terms (probably we could delete this attribute and
# temporary inlitialized in each child class.
# Then in each child class the specific terms arre added to the propagator delegate, and the propagation is done
#cycling through all the terms in propagate_one_step (for funct in propagator: func(..))

class Propagator():
    def __init__(self):
        self.propagator_terms = PropagatorTerms()
        self.propagator = []

    def init_propagaror_terms(self, dt, molecule, pcm):
        self.propagator_terms.set_attributes(dt, molecule, pcm)
        self.propagator_terms.init_terms_dictionary()

    def add_term_to_propagator(self, term_name):
            self.propagator.append(self.propagator_terms.dict_terms[term_name])

    def clean_propagator(self):
        self.propagator = []

    def set_propagator(self, dt, molecule, env):
        return

    def propagate_one_step(self, *args):

        return

    def propagate_n_step(self, *args):
        return










