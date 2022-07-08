from dataclasses import dataclass
from matplotlib import pyplot as plt
from numpy import ndarray
import numpy as np
from typing import Optional

@dataclass
class EnzymeKinetics:
    time: list
    substrate: Optional[list] = None
    product: Optional[list] = None
    inhibitor: Optional[list] = None
    init_substrate: Optional[list] = None

    def __post_init__(self):
        # TODO consisgtency check shape of input data

        self._is_substrate = self.check_is_substrate()
        self._multiple_concentrations = self._check_multiple_concentrations()
        if self.substrate is None:
            self.substrate = self.calculate_substrate()
        


    def check_is_substrate(self) -> bool:
        if self.substrate is not None:
            _is_substrate = True
        else:
            _is_substrate = False

        return _is_substrate



    def _check_multiple_concentrations(self) -> bool:
        """Checks if data contains one or multiple concentration array based on the shape of the array"""
        
        if self.substrate is not None and len(self.substrate.shape) == 2 or self.product is not None and len(self.product.shape) == 2:
            return True
        else:
            return False



    def calculate_substrate(self) -> ndarray:
        """If substrate data is not provided substrate data is calculated, assuming conservation of mass"""

        if self.substrate is None and self.product is not None:
            substrate = np.zeros(self.product.shape)
            if not self._multiple_concentrations:
                substrate = np.array([self.init_substrate - product for product in self.product])
            else:
                for i, row in enumerate(self.product):
                    substrate[i] = [self.init_substrate[i] - product for product in row]
                    #TODO: catch error if no init_substrate is provided
            
            return substrate

        else:
            raise Exception("Data must be provided eighter for substrate or product")


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from pyenzymekinetics.helper.load_utitlity import ek
