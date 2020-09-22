from OC.OCGeneticIterator import OCGeneticIterator
from OC.OCRabitzIterator import OCRabitzIterator
from OC.OCSimplePropagatorIterator import SimplePropagationIterator
from read_and_set.read.ReadGeneticConf import ReadGeneticConf
from read_and_set.read.ReadNoConf import ReadNoConf
from read_and_set.read.ReadRabitzConf import ReadRabitzConf
from read_and_set.set.SetGeneticOCInput import SetGeneticOCInput
from read_and_set.set.SetNoConf import SetNoConf
from read_and_set.set.SetRabitzInput import SetRabitzOCInput

OCAlgorithmDict = { "none": SimplePropagationIterator,
                    "genetic": OCGeneticIterator,
                    "rabitzi": OCRabitzIterator,
                    "rabitzii": OCRabitzIterator
}
OCAlgorithmConfig = {"none": ReadNoConf,
                    "genetic": ReadGeneticConf,
                    "rabitzi": ReadRabitzConf,
                    "rabitzii": ReadRabitzConf
                     }
OCAlgorithmSet = {"none": SetNoConf,
                    "genetic": SetGeneticOCInput,
                    "rabitzi": SetRabitzOCInput,
                    "rabitzii": SetRabitzOCInput
}
OCConfigFileDefaultNames = {"none": None,
                          "genetic": "genetic.conf",
                          "rabitzi": "rabitz.conf",
                          "rabitzii": "rabitz.conf"
                            }