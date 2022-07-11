from dataclasses import dataclass
from typing import Dict, Optional, Callable
from lmfit import Model, Parameters
from numpy import exp, ndarray


@dataclass
class CalibrationModel():
    name: str
    equation: Callable
    parameters: Optional[Dict[str, float]] = None
    lmfit_model: Optional[Model] = None
    lmfit_params: Optional[Parameters] = None
    result: Optional[Parameters] = None

    def __post_init__(self):
        self.lmfit_model = Model(self.equation)

    def set_lmfit_parameters(self) -> Parameters:
        return self.lmfit_model.make_params(**self.parameters)


# Model equations
# If change here, change also in utility.py

def linear1(x, a) -> float:
    return a*x


def quadratic(x, a, b) -> float:
    return a*x**2 + b*x


def poly3(x, a, b, c) -> float:
    return a*x**3 + b*x**2 + c*x


def poly_e(x, a, b) -> float:
    return a*exp(x/b)


def rational(x, a, b) -> float:
    return (a*x)/(b+x)
