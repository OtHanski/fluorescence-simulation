from SampleCell import SampleCell
from photon import photon
from photon import photonSave
import numpy as np

samples = 50
samp = SampleCell(np.linspace(0, 10, samples), np.ones(samples), samples=samples)

pos = np.array([0, 0, 5])

results = []

for i in range(10):
    try:
        phot = photon(sampCell=samp, position=pos)
        print(f"Simulating photon;\n{phot.pos}\n{phot.direction}\n")
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