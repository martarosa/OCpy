class FieldInput():
    def __init__(self):
        self.field_type = None
        # parameters used when field is created with a given shape
        self.dt = None
        self.nstep = None
        self.control_parameters = None
        self.fi = None
        self.omega = None
        self.sigma = None
        self.t0 = None

        self.omega_sys = None

        #if field is here assigned by values at each time and not created from parameters and shape
        self.field = None
        self.field_time_axis = None

        self.field_restart_name = None