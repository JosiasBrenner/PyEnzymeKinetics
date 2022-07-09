from logging.handlers import QueueHandler
from pyenzymekinetics.calibration.calibrationmodel import CalibrationModel, linear1, quadratic, poly3, poly_e, rational

from typing import Optional
from dataclasses import dataclass

from lmfit import minimize, Parameters, Parameter, report_fit, Model
from scipy.optimize import curve_fit
from numpy import ndarray




@dataclass
class StandardCurve():
    concentration: ndarray
    absorption: ndarray
    models: Optional[dict] = None

    def __post_init__(self):
        self.models = self.initialize_models()
        self.fit_models()


    def initialize_models(self) -> list:
        linear_model = CalibrationModel(
            name = "Linear",
            equation = linear1,
            parameters = {"a":0}
        )

        quadratic_model = CalibrationModel(
            name = "Quadratic",
            equation = quadratic,
            parameters = {"a":0, "b":0}
        )

        poly3_model = CalibrationModel(
            name = "3rd degree polynominal",
            equation = poly3,
            parameters = {"a":0, "b":0, "c":0}
        )

        polye_model = CalibrationModel(
            name = "Exponential polynominal",
            equation = poly_e,
            parameters = {"a":0, "b":0}
        )

        rational_model = CalibrationModel(
            name = "Rational",
            equation = rational,
            parameters = {"a":0, "b":0}
        )

        models = {linear_model.name: linear_model,
                  quadratic_model.name: quadratic_model,
                  poly3_model.name: poly3_model,
                  polye_model.name: polye_model,
                  rational_model.name: rational_model
                  }

        return models

    def fit_models(self):
        for model in self.models.values():
            for parameter in model.parameters.values():
                model.parameters = curve_fit(f=model.equation, xdata=self.concentration, ydata=self.absorption)[0]
            model.lmfit_params = CalibrationModel.set_lmfit_parameters(model)




if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from pyenzymekinetics.parameterestimator.helper.load_utitlity import calibration_conc, calibration_abso
    print(calibration_abso)
    print(calibration_conc)
    obj = StandardCurve(calibration_abso, calibration_conc)
    print("hi")
    
