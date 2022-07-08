import numpy as np 
from pyenzymekinetics import EnzymeKinetics

concentration_data = np.fromfile("data/concentration")
concentration_data = np.reshape(concentration_data,(7,21))/1000
time_data = np.fromfile("data/time")
init_substrate = np.array([1,2.5,5,7.5,10,20,30])
ek = EnzymeKinetics(time_data, product=concentration_data, init_substrate=init_substrate)
