from lmfit import Parameters
from typing import Dict, List
from pyenzymekinetics.parameterestimator.enzymekinetics import EnzymeKinetics

model_params_dict: Dict[str, List[str]] = {
    "irrev MM": ["kcat", "Km"],
    "irrev MM with enzyme inactivation": ["kcat", "Km", "ki"]
}


class KineticModel():
    def __init__(self) -> None:

    def set_params(self):
        kcat = self.get_initial_Km()/self.enzyme_conc
        km = self.get_initial_Km()

        params = Parameters()
        params.add('k_cat', value=kcat, min=kcat/100, max=kcat*100)
        params.add('Km', value=km, min=km/100,
                   max=np.max(self.substrate_conc)*100)

    def menten_irreversible(self, w, t, params):
        c_S, c_E = w

        k_cat = params['k_cat'].value
        K_m = params['Km'].value

        dc_S = k_cat * c_E * (c_S) / (K_m+c_S)
        dc_E = 0

        return (dc_S, dc_E)

    def set_startvector(self):
        return (self.init_substrate, self.enzyme)

    def irrev_MM(self):
        pass


if __name__ == "__main__":
    print("hi")
