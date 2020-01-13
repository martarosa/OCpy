import  f90nml
import sys
from read import NameListSections as sec

class ReadTDPLASNamelist():
    def __init__(self):
        self.TDPprop = sec.SectionTDPlasPropagate()
        self.TDPsur = sec.SectionTDPlasSurface()
        self.TDPeps = sec.SectionTDPlasEps()
        self.TDPmedium = sec.SectionTDPlasMedium()


    def read_file(self, folder, namefile):
        default_dict = {'propagate': {

                                    'interaction_stride': 1,
                                    'interaction_type':'pcm',
                                    'interaction_init':'non-scf',
                                    'propagation_type':'ief',
                                    'local_field':'loc',
                                    'out_level': 'low',
                                    'debug_type': 'non',
                                    'test_type': 'non',
                                    'scf_mix_coeff':0.2,
                                    'scf_max_cycles':600,
                                    'scf_threshold':10,
                                    'ntst':150,
                                    'print_lf_matrix':'non',
                                    'medium_relax': "non",},
                              'surface': {
                                    'input_surface': 'fil',
                                    'spheres_number':0,
                                    'sphere_radius': 0.0,
                                    'sphere_position_x':0.0,
                                    'sphere_position_y' : 0.0,
                                    'sphere_position_z':0.0,
                                    'spheroids_number':0,
                                    'spheroid_axis_x':0.0,
                                    'spheroid_axis_y': 0,
                                    'spheroid_axis_z': 0.0,
                                    'spheroid_position_x': 0.0,
                                    'spheroid_position_y': 0.0,
                                    'spheroid_position_z': 0.0,
                                    'spheroid_radius':0.0},
                              'medium': {
                                     'medium_init': 'fro',
                                     'medium_type': 'nan',
                                     'medium_pol': 'chr',
                                     'bem_type': 'diag',
                                     'bem_read_write': 'rea'},
                              'eps_function': {
                                     'epsilon_omega':'drl',
                                     'eps_0':1000.0,
                                     'eps_d':1.0,
                                     'eps_a':0.110224,
                                     'eps_gm':0.000757576,
                                     'eps_w0':0.0,
                                     'f_vel':0.0,
                                     'tau_deb':1000.0}
                              }

        tdplas_input = f90nml.read(folder + namefile)
        self.check_input_sections(tdplas_input)
        self.set_dictionaries(tdplas_input, default_dict)

    def set_dictionaries(self, tdplas_input, default_dict):
        self.TDPprop.init(tdplas_input, default_dict, "propagate")
        self.TDPprop.check()
        self.TDPmedium.init(tdplas_input, default_dict, "medium")
        self.TDPmedium.check()
        self.TDPeps.init(tdplas_input, default_dict, "eps_function")
        self.TDPeps.check()
        if tdplas_input["surface"]:
            print("Surface namelist is not used and takes only default values. \n"
                  +"Any value of this namelist is discharged\n")


    def check_input_sections(self, tdplas_input):
        sections = ["propagate" , "medium", "surface", "eps_function"]
        if not all(elem in tdplas_input.keys() for elem in sections):
            sys.exit("Error. Sections names are wrong in input file")



