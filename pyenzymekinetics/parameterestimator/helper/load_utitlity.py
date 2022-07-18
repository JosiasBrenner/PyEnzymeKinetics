import numpy as np
#from pyenzymekinetics.parameterestimator import EnzymeKinetics

absorbance_measured = np.fromfile("data/concentration")
absorbance_measured = np.reshape(absorbance_measured, (7, 21))
calibration_abso = np.fromfile("data/calibration_abso")
calibration_conc = np.fromfile("data/calibration_conc")
time = np.linspace(0, 20, 21)

time_data = np.fromfile("data/time")
init_substrate = np.array([1, 2.5, 5, 7.5, 10, 20, 30])

product_chantal = np.fromfile("data/example_chantal/product")
time_chantal = np.fromfile("data/example_chantal/time")
init_sub_chantal = np.fromfile("data/example_chantal/init_sub")

chantal = (product_chantal, time_chantal, init_sub_chantal)
