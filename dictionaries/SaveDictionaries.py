from save.SaveOCGenetic import SaveOCGenetic
from save.SaveOCRabitz import SaveOCRabitz
from save.SaveSimplePropagation import SaveSimplePropagation
from save.SaveOCScipyOptimize import SaveOCScipyOptimize

SaveDict = {
    "genetic": SaveOCGenetic,
    "rabitzi": SaveOCRabitz,
    "rabitzii":SaveOCRabitz ,
    "none": SaveSimplePropagation,
    "nelder-mead": SaveOCScipyOptimize,
    "cobyla": SaveOCScipyOptimize,
    "bfgs": SaveOCScipyOptimize,
    "cg": SaveOCScipyOptimize
}



