from medium.DinamicMedium import DinamicMedium
from medium.FrozenSolventMedium import FrozenSolventMedium
from medium.VacMedium import VacMedium

MediumDict = {
    "vac": VacMedium,
    "sol": FrozenSolventMedium,
    "nanop": DinamicMedium
}