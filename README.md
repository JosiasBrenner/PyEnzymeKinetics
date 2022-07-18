# PyEnzymeKinetics

Software-package for easy enzyme kinetic parameter estimation.
Supports EnzymeML format (not yet)
Allows comprehensive data analysis from experimental raw data to kinetic parameters

## Calibration

- linear and non-linear calibration for calculating concentrations from analytic signal (spectroscopy / HPLC)
- Selection of best fit model based on akaike-criterion

## Michaelis Menten kinetics

- Fitting of enzyme assay data Michaelis Menten equations for:
  - irreversible MM
  - enzyme inactivation
  - inhibition

### Dependencies

- standard stuff

### Installation

```
pip install PyEnzymeKinetics
```

### Demo
By providing a list of concetrations and the respective measured absorptions, the data is fitted against multiple equation by running:
```
from pyenzymekinetics import CalibrationModel
calibration = CalibrationModel(concentration_array, absorption_array)
```
The resulting standard curve model can be visualiuzed by calling:
```
calibration.visualize()
```

## Help

Any advise for common problems or issues.

```
command to run if program contains helper info
```

## Authors

Contributors names and contact info

Max Haeussler
max.haeussler@t-online.de


## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments

ex
