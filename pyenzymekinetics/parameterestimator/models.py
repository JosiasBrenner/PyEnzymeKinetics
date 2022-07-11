from lmfit import Parameters
from typing import Dict


class KineticModel:
    def __init__(self, is_substrate: bool) -> None:
        self.is_substrate = is_substrate
        self.params = None,
        self.equation = None
        pass

    def irrev_MM(self):
        pass


hallo = KineticModel(True)

hmm: Dict[KineticModel, str] = {
    hallo: "max"
}

print(hmm)
