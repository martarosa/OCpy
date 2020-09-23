from read_and_set.read.field_restart.ReadGeneralFieldRestart import ReadGeneralFieldRestart
from read_and_set.read.field_restart.ReadGeneticFieldRestart import ReadGeneticFieldRestart



FieldRestartDict = {
    "eulero_2order": ReadGeneralFieldRestart,
    "eulero_1order": ReadGeneralFieldRestart,
    "rabitzi": ReadGeneralFieldRestart,
    "rabitzii": ReadGeneralFieldRestart,
    "genetic": ReadGeneticFieldRestart
}
