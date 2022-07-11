from pyenzymekinetics.calibration.standardcurve import StandardCurve

from typing import Dict, Callable
from numpy import ndarray, zeros, exp
from scipy.optimize import fsolve


# Model equations
# If change here, change also in calibrationmodel.py
def root_linear1(x: float, params: Dict[str, float]) -> float:
    a, absorption = params.values()
    return a*x - absorption


def root_quadratic(x: float, params: Dict[str, float]) -> float:
    a, b, absorption = params.values()
    return a*x**2 + b*x - absorption


def root_poly3(x: float, params: Dict[str, float]) -> float:
    a, b, c, absorption = params.values()
    return a*x**3 + b*x**2 + c*x - absorption


def root_poly_e(x: float, params: Dict[str, float]) -> float:
    a, b, absorption = params.values()
    return a*exp(x/b) - absorption


def root_rational(x: float, params: Dict[str, float]) -> float:
    a, b, absorption = params.values()
    return (a*x)/(b+x) - absorption

# Mapper for equations from CalibrationModel
equation_dict: Dict[str, Callable] = {
    "Linear": root_linear1,
    "Quadratic": root_quadratic,
    "3rd polynominal": root_poly3,
    "Exponential": root_poly_e,
    "Rational": root_rational
    }

def to_concentration(standard_curve: StandardCurve, data: ndarray, standard_curve_model="") -> ndarray:
    """Converts absobrance / peak area into concentration based on StandardCurve.

    Args:
        standard_curve (StandardCurve): StandardCurve object, containing model for 
        concentration calculation.

        data (ndarray): Measured peak area / absorbance

        standard_curve_model (str, optional): Optionally, a differing model from the best model 
        from StandardCurve can be provided by passing the name of the respective model. Defaults to "".

    Returns:
        ndarray: Calculated concentrations
    """
    
    result = zeros(data.shape)

    best_model = next(iter(standard_curve.result_dict))
    if len(standard_curve_model) == 0:
        model = best_model
    else:
        model = standard_curve_model

    equation: Callable = equation_dict[standard_curve.models[model].name]
    params: dict = standard_curve.models[model].result.params.valuesdict()

    measurments, values = data.shape
    for measurement in range(measurments):
        for value in range(values):
            params["absorption"] = data[measurement, value]
            result[measurement, value] = fsolve(equation, 0, params)

    return result

if __name__ == "__main__":
    from pyenzymekinetics.parameterestimator.helper.load_utitlity import absorbance_measured, calibration_abso, calibration_conc
    print(absorbance_measured.shape)
    calibration = StandardCurve(calibration_conc, calibration_abso, "mM")

    print(to_concentration(calibration, absorbance_measured))
    