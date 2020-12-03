from OC.OCGeneticIterator import OCGeneticIterator
from OC.OCRabitzIterator import OCRabitzIterator
from OC.OCSimplePropagatorIterator import SimplePropagationIterator
from OC.OCScipyIterator import OCScipyOptimizeIterator
from ObjectiveFunction.ObjectiveFunctionSpecifics import ObjectiveFunctionGroundState, ObjectiveFunctionOptical
from read_and_set.read.conf.ReadGeneticConf import ReadGeneticConf
from read_and_set.read.conf.ReadNoConf import ReadNoConf
from read_and_set.read.conf.ReadRabitzConf import ReadRabitzConf
from read_and_set.set.SetGeneticConf import SetGeneticOCInput
from read_and_set.set.SetNoConf import SetNoConf
from read_and_set.set.SetRabitzConf import SetRabitzOCInput

OCObjectiveFunction = {"optical_excitation": ObjectiveFunctionOptical,
                 "ground_state": ObjectiveFunctionGroundState}


OCAlgorithmDict = { "none": SimplePropagationIterator,
                    "genetic": OCGeneticIterator,
                    "rabitzi": OCRabitzIterator,
                    "rabitzii": OCRabitzIterator,
                    "nelder-mead": OCScipyOptimizeIterator,
                    "bfgs": OCScipyOptimizeIterator,
                    "cg": OCScipyOptimizeIterator,
                    "cobyla": OCScipyOptimizeIterator
}


OCAlgorithmConfig = {"none": ReadNoConf,
                     "nelder-mead": ReadNoConf,
                     "bfgs": ReadNoConf,
                     "cg": ReadNoConf,
                    "genetic": ReadGeneticConf,
                    "rabitzi": ReadRabitzConf,
                    "rabitzii": ReadRabitzConf,
                    "cobyla": ReadNoConf
                     }
OCAlgorithmSet = {  "none": SetNoConf,
                  "nelder-mead": SetNoConf,
                  "bfgs": SetNoConf,
                  "cg": SetNoConf,
                    "genetic": SetGeneticOCInput,
                    "rabitzi": SetRabitzOCInput,
                    "rabitzii": SetRabitzOCInput,
                    "cobyla": SetNoConf
}


OCConfigFileDefaultNames = {"none": None,
                            "nelder-mead": None,
                            "cg": None,
                            "bfgs": None,
                            "cobyla": None,
                            "genetic": "genetic.conf",
                            "rabitzi": "rabitz.conf",
                            "rabitzii": "rabitz.conf"
                            }