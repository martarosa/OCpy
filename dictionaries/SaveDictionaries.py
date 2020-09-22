from save.SaveOCGenetic import SaveOCGenetic
from save.SaveOCRabitz import SaveOCRabitz
from save.SaveSimplePropagation import SaveSimplePropagation

SaveDict = {
    "genetic": SaveOCGenetic,
    "rabitzi": SaveOCRabitz,
    "rabitzii":SaveOCRabitz ,
    "none": SaveSimplePropagation,
}



