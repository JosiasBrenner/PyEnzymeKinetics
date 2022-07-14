from pyenzyme import EnzymeMLDocument
import numpy as np


def get_v(substrate_conc, time):
    v_all = 0.0*substrate_conc[:]  # initialize velocity vector
    if len(substrate_conc.shape) > 1:
        for i in range(substrate_conc.shape[0]):

            prev_value = substrate_conc[i, 0]
            prev_time = 0.0

            for j in range(substrate_conc.shape[1]):

                if time[j] == 0:
                    delta = prev_value - substrate_conc[i, j]
                else:
                    delta = abs(
                        (prev_value - substrate_conc[i, j])/(time[j]-prev_time))

                v_all[i, j] = delta
                prev_value = substrate_conc[i, j]
                prev_time = time[j]

        v = np.max(v_all, axis=0)

    else:

        prev_value = substrate_conc[0]
        prev_time = 0.0

        for j in range(substrate_conc.shape[0]):

            if time[j] == 0:
                delta = prev_value - substrate_conc[j]
            else:
                delta = abs(
                    (prev_value - substrate_conc[j])/(time[j]-prev_time))

            v_all[j] = delta
            prev_value = substrate_conc[j]
            prev_time = time[j]

        v = v_all

    return v


def writeConcentrationData(enzmldoc: EnzymeMLDocument, species_ID: str, data: np.array, data_unit=""):
    """overwrites measurements of EnzymeML document with array.

    Args:
        enzmldoc (EnzymeMLDocument): document, which should be overwritten
        species_ID (str): species of document (e.g 's0')
        data (np.array): array, shaped according to the data within the EnzymeML document
        data_unit (str, optional): unit of the data provided in the array. Defaults to "".
    """
    for i, measurement in enumerate(enzmldoc.measurement_dict.keys()):
        enzmldoc.getMeasurement(measurement).getReactant(
            species_ID).replicates[0].data = list(data[i])
        enzmldoc.getMeasurement(measurement).getReactant(
            species_ID).replicates[0].data_type = "conc"
        enzmldoc.getMeasurement(measurement).getReactant(
            species_ID).replicates[0].is_calculated = True
        if len(data_unit) != 0:
            enzmldoc.getMeasurement(measurement).getReactant(
                species_ID).replicates[0].data_unit = data_unit
