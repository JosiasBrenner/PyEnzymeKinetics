from pyenzymekinetics.utility.initial_parameters import get_v, get_initial_vmax, get_initial_K

from matplotlib import pyplot as plt
from numpy import ndarray, array, zeros, max


class EnzymeKinetics():

    def __init__(self,
                 time: ndarray,
                 enzyme: ndarray,
                 substrate: ndarray = None,
                 product: ndarray = None,
                 init_substrate: ndarray or float = None,
                 inhibitor: ndarray = None
                 ):
        self.time = time
        self.enzyme = enzyme
        self.substrate = substrate
        self.product = product
        self.init_substrate = init_substrate
        self.inhibitor = inhibitor

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
            substrate = zeros(self.product.shape)
            if not self._multiple_concentrations:
                substrate = array(
                    [self.init_substrate - product for product in self.product])
            else:
                for i, row in enumerate(self.product):
                    substrate[i] = [self.init_substrate[i] -
                                    product for product in row]
                    # TODO: catch error if no init_substrate is provided

            return substrate

        else:
            raise Exception(
                "Data must be provided eighter for substrate or product")


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from pyenzymekinetics.parameterestimator.helper.load_utitlity import *
    from pyenzymekinetics.calibrator.standardcurve import StandardCurve
    from pyenzymekinetics.calibrator.utility import to_concentration

    # Calibrate
    standardcurve = StandardCurve(calibration_conc, calibration_abso)
    # standardcurve.visualize_fit()

    # Convert concentration in absorbance data
    conc = to_concentration(standardcurve, absorbance_measured)

    kinetics = EnzymeKinetics(
        time, enzyme=0.8, product=conc, init_substrate=init_substrate)

    print(kinetics._get_initial_vmax())
    print("hi")
