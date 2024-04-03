from SampleCell import SampleCell
from photon import photon
from photon import photonSave
import numpy as np


samp = SampleCell(np.linspace(0, 10, 100), np.ones(100), samples=100)

results = []

for i in range(10):
    phot = photon(samp)
    result = phot.simulate()
    results.append(result)



#ERROR
#line 126, in hit_wall
#    p = [self.specular_probability[wavelength][z_index],
#         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^
#IndexError: index 99 is out of bounds for axis 0 with size 99
    


photonSave("./testData.csv", results)