class LogHeader():
    def __init__(self):
        self.header = None
        self.field_header = None
        self.oc_header = None
        self.oc_conf_header =  None

    def init_log_header(self, log_input):
        self.init_oc_header(log_input)
        if log_input.restart == "norestart_found":
            restart = ("#WARNING: restart asked but restart file not found. Starting from default field" + "\n")
        elif log_input.restart == "true":
            restart = ("#Restarted from field \n")
            if log_input.oc_algorithm == "rabitzi" or log_input.oc_algorithm == "rabitzii":
                log_input.field_type = "restart_rabitz"
            elif log_input.oc_algorithm == "genetic":
                log_input.field_type = "restart_genetic"
        elif log_input.restart == "only_bkp_found":
            restart = ("#WARNING: restart asked but restart file not found. Restarting from bkp field \n")
            if log_input.oc_algorithm == "rabitzi" or log_input.oc_algorithm == "rabitzii":
                log_input.field_type = "restart_rabitz"
            elif log_input.oc_algorithm == "genetic":
                log_input.field_type = "restart_genetic"
        else:
            restart = "\n"
        self.init_field_header(log_input)
        self.init_conf_header(log_input)
        print(self.oc_conf_header)
        self.header = ( "#calculation: " + log_input.oc_algorithm + "\n" +
                        self.oc_header +
                        self.oc_conf_header +
                        "#dt: " + log_input.dt + " medium: " + log_input.medium + "\n" +
                        restart +
                        self.field_header)



    def init_oc_header(self, log_input):
        if (log_input.oc_algorithm != "none"):
            self.oc_header = ("#target state: " + log_input.target_state + "\n" +
                              "#alpha: " + log_input.alpha + "\n\n"
                              )
        elif (log_input.oc_problem == "ground_state"):
            self.oc_header = ("#Ground state computation in the CISD determinant space from HF reference wavefunction")
        else:
            self.oc_header="\n"


    def init_conf_header(self, conf_input):
        if (conf_input.oc_algorithm not in  ["none", "nelder-mead", "bfgs", "cg", "cobyla"] ):
            self.oc_conf_header = ("#oc algorithm input: " + conf_input.string_conf_file + "\n\n")
        else:
            self.oc_conf_header = "\n"


    def init_field_header(self, log_input):
        field = {
            'const': lambda: self.const_pulse(log_input),
            'pip': lambda: self.pip_pulse(log_input),
            'sin': lambda: self.sin_pulse(log_input),
            'gau': lambda: self.gau_pulse(log_input),
            'sum': lambda: self.sum_pulse(log_input),
            'genetic': lambda: self.genetic_pulse(log_input),
            'gaussian_sum': lambda: self.gaussian_sum(log_input),
            'free_harmonics': lambda: self.free_harmonics(log_input),
            # only internal values
            'restart_rabitz': lambda: self.restart_rabitz(log_input),
            'restart_genetic': lambda: self.genetic_pulse(log_input)
        }
        return field.get(log_input.field_type, lambda: "Inexistent field type")()




    def restart_rabitz(self, log_input):
        self.field_header = ""
        pass
    
    def free_harmonics(self, log_input):
        self.field_header = ("#field parameters: \n#field: " + log_input.field_type + "\n"
                             "#number of control parameters: " + log_input.num_control_parameters + "\n"
                             "#omega: initial guess field is no perturbation active or initialized at random \n\n")     
    
    def gaussian_sum(self, log_input):
        self.field_header = ("#field parameters: \n#field: " + log_input.field_type + "\n"
                             "#number of control parameters: " + log_input.num_control_parameters + "\n" 
                             "#omega: initial guess field is no perturbation active or initialized at random \n\n")
                             


    def const_pulse(self, log_input):
        self.field_header = ("#field parameters: \n#field: " + log_input.field_type + "\n" 
                             "#fi: " + log_input.fi +"\n\n"
                             )

    def pip_pulse(self, log_input):
        self.field_header = ("#field parameters: \n#field: " + log_input.field_type +"\n"
                             "#fi: " + log_input.fi + "\n" +
                             "#omega: " + log_input.omega + "\n"
                             "#sigma: " + log_input.sigma +" t0: " + log_input.t0 + "\n\n")


    def sin_pulse(self, log_input):
        self.field_header = ("#field parameters: \n#field: " + log_input.field_type + "\n"
                             "#fi: " + log_input.fi + "\n" +
                             "#omega: " + log_input.omega + "\n"
                             "#sigma: " + log_input.sigma + "\n\n")

    def gau_pulse(self, log_input):
        self.field_header = ("#field parameters: \n#field: " + log_input.field_type + "\n"
                             "#fi: " + log_input.fi + "\n" +
                             "#sigma: " + log_input.sigma + " t0: " + log_input.t0 + "\n\n")

    def sum_pulse(self, log_input):
        self.field_header = ("#field parameters: \n#field: " + log_input.field_type + "\n"
                             "#fi: " + log_input.fi + "\n" +
                             "#omega: " + log_input.omega + "\n\n")

    def genetic_pulse(self, log_input):
        self.field_header = ("#field parameters: \n#field: " + log_input.field_type + "\n" 
                             "#omega: fourier frequencies \n\n")



