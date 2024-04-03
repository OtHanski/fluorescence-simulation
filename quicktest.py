import SampleCell as sc
import photon as ph
import numpy as np

samp = sc.SampleCell(np.linspace(0, 10, 100), np.ones(100), samples=100)

results = []

for i in range(10):
    phot = ph.photon(samp)
    result = phot.simulate()
    results.append(result)

print(results)


#ERROR
#line 126, in hit_wall
#    p = [self.specular_probability[wavelength][z_index],
#         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^
#IndexError: index 99 is out of bounds for axis 0 with size 99