from lmfit import Parameters
from typing import Dict, List, Callable, Tuple
from pyenzymekinetics.utility.initial_parameters import get_v, get_initial_Km, get_initial_vmax
from numpy import max

model_params_dict: Dict[str, List[str]] = {
    "irrev MM": ["kcat", "Km"],
    "irrev MM with enzyme inactivation": ["kcat", "Km", "ki"]
}


class KineticModel():
    def __init__(self,
                 name: str,
                 model: Callable,
                 params: str,
                 kcat_initial: float,
                 Km_initial: float
                 ) -> None:

        self.name = name
        self.model = model
        self.params = params
        self.kcat_initial = kcat_initial
        self.Km_initial = Km_initial
        self.parameters = self.set_params(*self.params)
        self.w0 = None
        self.result = None

    def set_params(self, *args):

        params = Parameters()
        params.add('k_cat', value=self.kcat_initial,
                   min=self.kcat_initial/100, max=self.kcat_initial*100)
        params.add('Km', value=self.Km_initial, min=self.Km_initial/100,
                   max=max(self.Km_initial)*1000)

        if self.params == "ki":
            params.add(self.params, value=0.01, min=0.0001, max=0.9999)

        return params


def menten_irreversible(w0: tuple, t, params):
    c_S, c_E = w0

    k_cat = params['k_cat'].value
    K_m = params['Km'].value

    dc_S = -k_cat * c_E * (c_S) / (K_m+c_S)
    dc_E = 0

    return (dc_S, dc_E)


def menten_irreversible_enzyme_inact(w0: tuple, t, params) -> tuple:
    c_S, c_E = w0

    k_cat = params['k_cat'].value
    K_m = params['Km'].value
    k_i = params["ki"].value

    dc_S = -k_cat * c_E * (c_S) / (K_m+c_S)
    dc_E = -k_i * c_E

    return (dc_S, dc_E)


if __name__ == "__main__":
    print("hi")
