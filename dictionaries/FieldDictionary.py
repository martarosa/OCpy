from read_and_set.read.field_restart.ReadGeneralFieldRestart import ReadGeneralFieldRestart
from read_and_set.read.field_restart.ReadGeneticFieldRestart import ReadGeneticFieldRestart



FieldRestartDict = {
    "none": ReadGeneralFieldRestart,
    "rabitzi": ReadGeneralFieldRestart,
    "rabitzii": ReadGeneralFieldRestart,
    "genetic": ReadGeneticFieldRestart
}
