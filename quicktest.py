from SampleCell import SampleCell
from photon import photon
from photon import photonSave
import numpy as np

specrefl = {"121.567E-9": np.zeros(1000-1)+1, "450E-9": np.zeros(1000-1)+0.98}
diffrefl = {"121.567E-9": np.zeros(1000-1), "450E-9": np.zeros(1000-1)}
absprob = {"121.567E-9": np.zeros(1000-1), "450E-9": np.zeros(1000-1)+0.02}
WLconversion = {"121.567E-9": np.zeros(1000-1), "450E-9": np.zeros(1000-1)}

samp = SampleCell(np.linspace(0, 10, 1000), np.ones(1000), samples=1000, specrefl=specrefl, diffrefl=diffrefl, absprob=absprob, WLconversion=WLconversion)

pos = np.array([0, 0, 5])

results = []
negdirAfter = 0
negdirBefore = 0
totalPhot = 0
totalPhot2 = 0
dirrrr: np.ndarray = np.array([0.30960542, -0.09976733, -0.94561671])
negdirs1 = []
negdirs2 = []
for i in range(100):
    phot = photon(sampCell=samp, position=pos)
    negdirs1.append(phot.getDir())
    try:
        totalPhot += 1
        if phot.getDir()[2] < 0:
            negdirBefore += 1
        result = phot.simulate()
        results.append(result)
        if phot.getDir()[2] < 0:
            negdirAfter += 1
        totalPhot2 += 1
        negdirs2.append(phot.getDir())
    except:
        negdirs2.append(420.69)
        continue
print(negdirs1[0], "\n", negdirs2[0])

print("z neg before:", negdirBefore, "\nz neg after:", negdirAfter, "\nPercentage died:", 100 - negdirAfter/negdirBefore*100)
print("total:", totalPhot, totalPhot2)

#ERROR
#line 126, in hit_wall
#    p = [self.specular_probability[wavelength][z_index],
#         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^^^^^^^^
#IndexError: index 99 is out of bounds for axis 0 with size 99
    


#photonSave("./data/testData11_2.csv", results)