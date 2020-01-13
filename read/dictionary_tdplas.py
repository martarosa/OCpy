class TDPlas_legend_elements():
    def __init__(self, user_key, internal_key, user_values, internal_values):
        self.user_input_key = user_key
        self.internal_key = internal_key
        self.user_input_values = user_values
        self.internal_values = internal_values



class TDPlas_legend():
    def __init__(self):
        self.legend = []


    def init(self):
        legend=[
            ["interaction_stride", "n_q", [], []],
            ["interaction_init", "Finit_int", ["scf_es", "non-scf"], ["sce","nsc"]],
            ["interaction_type", "Fint", ["ons", "pcm"], ["ons", "pcm"]],
            ["propagation_type", "Fprop", ["ief", "ied", "ons","dip"], ["chr-ief","chr-ied","chr-ons","dip"]],
            ["scf_mix_coeff", "mix_coef", [], []],
            ["scf_threshold", "thrshld", [],[]],
            ["scf_max_cycles", "ncymax", [], []],
            ["local_field", "Floc", ["loc"," "], ["loc", " "]],
            ["out_level", "Fwrite", ["high", "low"], ["high", "low"]],
            ["debug_type", "Fdeb", ["equ", "vmu", "off", "non"], ["equ", "vmu", "off", "non"]],
            ["test_type", "Ftest", ["n-r", "n-l", "s-r", "s-l", "qmt","non"],["n-r", "n-l", "s-r", "s-l", "qmt","non"]],
            ["ntst", None,[],[]],
            ["medium_relax", "Fmdm_relax", ["rel", "non"], ["rel", "non"]],
            ["print_lf_matrix", None, [], []],
           ["input_surface", "Fsurf", ["fil", "gms","inv","non", "bui"], ["fil", "gms","inv","non", "bui"]],
            ['medium_init', "Finit_mdm", ['vac', "fro", "rea"], ['vac', "fro", "rea"]],
            ['medium_type',"Fmdm",['nan',"sol", "qso","qna"],["Cnan", "Csol", "Qsol", "Qnan"]],
            ['medium_pol', "Fmdm_pol", ['chr, "dip'], ['chr, "dip']],
            ['bem_type', "Fbem", ['diag', "stan"], ["diag", "stan"]],
            ['bem_read_write',"FinitBEM", ['rea', 'wri'], ['rea', 'wri']],
            ["epsilon_omega", "Feps", ["deb", "drl", "gen"], ["deb", "drl", "gen"]]

        ]


        for i in range(len(legend)):
            self.legend.append(TDPlas_legend_elements(legend[i][0], legend[i][1], legend[i][2], legend[i][3]))



