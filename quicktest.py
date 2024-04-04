from SampleCell import SampleCell
from photon import photon
from photon import photonSave
import numpy as np


samp = SampleCell(np.linspace(0, 10, 1000), np.ones(1000), samples=1000)

pos = np.array([0, 0, 5])

results = []

for i in range(10000):
    try:
        phot = photon(sampCell=samp, position=pos)
        result = phot.simulate()
        results.append(result)
    except:
        continue



#ERROR
#line 126, in hit_wall
#    p = [self.specular_probability[wavelength][z_index],
#         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^
#IndexError: index 99 is out of bounds for axis 0 with size 99
    


photonSave("./testData2.csv", results)