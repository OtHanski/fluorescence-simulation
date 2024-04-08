from SampleCell import SampleCell
from photon import photon
from dataAnalysis import saveData
import numpy as np

specrefl = {"121.567E-9": np.zeros(1000-1)+1, "450E-9": np.zeros(1000-1)+0.98}
diffrefl = {"121.567E-9": np.zeros(1000-1), "450E-9": np.zeros(1000-1)}
absprob = {"121.567E-9": np.zeros(1000-1), "450E-9": np.zeros(1000-1)+0.02}
WLconversion = {"121.567E-9": np.zeros(1000-1), "450E-9": np.zeros(1000-1)}

samp = SampleCell(np.linspace(0, 10, 1000), np.ones(1000), samples=1000)

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
        negdirs2.append("Good")
    except:
        negdirs2.append("Bad")
        continue
#print(negdirs1, "\n", negdirs2)

#print("z neg before:", negdirBefore, "\nz neg after:", negdirAfter, "\nPercentage died:", 100 - negdirAfter/negdirBefore*100)
#print("total:", totalPhot, totalPhot2)

good = []
bad = []
for i in range(len(negdirs1)):
    if negdirs2[i] == "Good":
        if negdirs1[i][2] < 0:
            good.append(negdirs1[i])
        else:
            continue
    else:
        bad.append(negdirs1[i])

#print("GOOD\n", good)
#print("BAD\n", bad)

    


#saveData("./data/testData3.csv", results)