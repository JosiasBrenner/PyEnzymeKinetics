from dataclasses import dataclass
from typing import Optional
from lmfit import Model, Parameters
from numpy import exp



@dataclass
class CalibrationModel():
    name: str
    equation: callable
    parameters: Optional[dict] = None
    lmfit_model: Optional[Model] = None
    lmfit_params: Optional[Parameters] = None


    def __post_init__(self):
        self.lmfit_model = Model(self.equation)
        #self.lmfit_params = self.set_lmfit_parameters()


    def set_lmfit_parameters(self) -> Parameters:
        return self.lmfit_model.make_params(**self.parameters)


# Model equations
def linear1(x, a):
    return a*x

def quadratic(x, a, b):
    return a*x**2 + b*x

def poly3(x, a, b, c):
    return a*x**3 + b*x**2 + c*x

def poly_e(x, a, b):
    return a*exp(x/b)

def rational(x, a, b):
    return (a*x)/(b+x)