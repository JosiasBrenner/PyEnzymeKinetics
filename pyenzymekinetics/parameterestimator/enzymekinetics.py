from opcode import haslocal
from pyenzymekinetics.utility.initial_parameters import get_v, get_initial_vmax, get_initial_Km
from pyenzymekinetics.parameterestimator.models import KineticModel, menten_irreversible, menten_irreversible_enzyme_inact, menten_irreversible_inhibition

from typing import Dict
from matplotlib import pyplot as plt
from numpy import ndarray, array, zeros, max
from scipy.integrate import odeint
from lmfit import minimize, report_fit


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

        self._is_substrate = self._check_is_substrate()
        self._multiple_concentrations = self._check_multiple_concentrations()
        if self.substrate is None:
            self.substrate = self.calculate_substrate()
        self._w0 = self._get_w0()

        self.models: Dict[str, KineticModel] = self.initialize_models()

    def _check_is_substrate(self) -> bool:
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

    def _get_w0(self):
        return (self.init_substrate, self.enzyme)

    def _get_kcat(self) -> float:
        return get_initial_vmax(self.substrate, self.time) / self.enzyme

    def initialize_models(self) -> Dict[str, KineticModel]:
        irreversible_MM = KineticModel(
            name="irreversible Michaelis Menten",
            params=(""),
            w0={"cS": self.init_substrate, "cE": self.enzyme, "cP": self.product},
            kcat_initial=self._get_kcat(),
            Km_initial=get_initial_Km(self.substrate, self.time),
            model=menten_irreversible
        )

        irrev_MM_enz_inact = KineticModel(
            name="irreversible Michaelis Menten with enzyme inactivation",
            params="ki",
            w0={"cS": self.init_substrate, "cE": self.enzyme, "cP": self.product},
            kcat_initial=self._get_kcat(),
            Km_initial=get_initial_Km(self.substrate, self.time),
            model=menten_irreversible_enzyme_inact
        )

        irrev_MM_prod_inhib = KineticModel(
            name="irreversible Michaelis Menten with competitive inhibition",
            params="kpi",
            w0={"cS": self.init_substrate, "cE": self.enzyme,
                "cP": self.product, "cI": self.inhibitor},
            kcat_initial=self._get_kcat(),
            Km_initial=get_initial_Km(self.substrate, self.time),
            model=menten_irreversible_inhibition
        )

        irrev_MM_prod_inhib_enz_inact = KineticModel(
            name="irreversible Michaelis Menten with competitive inhibition and enzyme inactivation",
            params=["kpi", "ki"],
            w0={"cS": self.init_substrate, "cE": self.enzyme,
                "cP": self.product, "cI": self.inhibitor},
            kcat_initial=self._get_kcat(),
            Km_initial=get_initial_Km(self.substrate, self.time),
            model=menten_irreversible_inhibition
        )

        kinetic_model_dict: Dict[str, KineticModel] = {
            irreversible_MM.name: irreversible_MM,
            irrev_MM_enz_inact.name: irrev_MM_enz_inact,
            irrev_MM_prod_inhib.name: irrev_MM_prod_inhib,
            irrev_MM_prod_inhib_enz_inact.name: irrev_MM_prod_inhib_enz_inact
        }

        return kinetic_model_dict

    def fit_models(self):
        for kineticmodel in self.models.values():

            def g(t, w0, params):
                '''
                Solution to the ODE w'(t)=f(t,w,p) with initial condition w(0)= w0 (= [S0])
                '''
                w = odeint(kineticmodel.model, w0, t, args=(params,))
                return w

            def residual(params, t, data):

                # get dimensions of data (here we fit against 4 measurments => ndata = 4)
                ndata, nt = data.shape
                resid = 0.0 * data[:]  # initialize the residual vector

                # compute residual per data set
                for i in range(ndata):

                    if len(kineticmodel.w0.keys()) == 3:
                        cS, cE, cP = kineticmodel.w0.values()
                        # TODO: fix initia product concentration
                        w0 = (cS[i], cE, 0)
                    else:
                        cS, cE, cP, cI = kineticmodel.w0.values()
                        w0 = w0 = (cS[i], cE, 0, cI)

                    model = g(t, w0, params)  # solve the ODE with sfb.

                    # get modeled product
                    model = model[:, 0]

                    # compute distance to measured data
                    resid[i, :] = data[i, :]-model

                return resid.flatten()

            kineticmodel.result = minimize(residual, kineticmodel.parameters, args=(
                self.time, self.substrate), method='leastsq', nan_policy='omit')

    def visualize(self, save_to_path="", format="svg"):
        for model in self.models.values():
            for i in range(self.substrate.shape[0]):
                ax = plt.scatter(x=self.time, y=self.substrate[i])

                # Integrate model
                s0 = self.init_conc[i]
                p0 = self.substrate_conc[i, 0]
                e0 = self.enzyme_conc

                #print(f"P: {p0}, E: {e0}, s:{s0}")

                w0 = (p0, e0, s0)

                data_fitted = self.g(self.time, w0, self.result.params)
                ax = plt.plot(self.time, data_fitted[:, 0], '-', linewidth=1)


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
    kinetics.fit_models()
    for model in kinetics.models.values():
        report_fit(model.result)

    print("hi")
