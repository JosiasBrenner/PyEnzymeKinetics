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
