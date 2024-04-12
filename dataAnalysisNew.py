import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from SampleCell import SampleCell
import FileHandler as fh

#========================SETUP=======================================================

fileName = ""
if not fileName:
    fileName = fh.ChooseSingleFile(initdir = "./data")
    print(fileName)

sampCellRadius = 5E-3
sampCellZ = 100E-3

top = 1 # Picks photons that exited at sample cell top
posDistributionPlot = 1
angleDistributionPlot = 1
xyPlanePlot = 1
wallHeatMapPlot = 1
# To be continued
numOfWallHitHist = 0
printInfo = 1

saveFigure = 0
posDistImName = "data/simDistTop.png"
angDistImName = "data/simAngTopSampleInTheMiddle.png"
xyPlaneImName = "data/simXYBotSampleInTheMiddle.png"
wallHeatMapImName = "data/wallHeatMapSampleInTheMiddle.png"

#===========================NOTES==========================================================================

# uv: specRefl 25%, absorp 25%, conversion 50%
# blue: specRefl 98%, absorp 2%
# 5 mm cell radius, 100 mm cell height
# 1 mm gas radius, 10 mm gas height

#========================FUNCTIONS======================================================================0===



# Read data from JSON
def readJsonData(fileName):
    stuff = {"cellSpecs": str, 
             "pos": np.ndarray, 
             "dir": np.ndarray, 
             "wallHits": np.ndarray, 
             "wavelength": np.ndarray,
             "event": np.ndarray,
             "angle": np.ndarray}
    if not fileName:
        fileName = fh.ChooseSingleFile(initdir = "./data")
    readdata = fh.ReadJson(fileName)
    pos, dir, wallHits, wavelength, event, angle = [], [], [], [], [], []
    for i in range(1,len(readdata["photons"])+1):
        pos.append(readdata["photons"][str(i)]["position"])
        dir.append(readdata["photons"][str(i)]["direction"])
        wallHits.append(readdata["photons"][str(i)]["bounces"])
        wavelength.append(readdata["photons"][str(i)]["wavelength"])
        event.append(readdata["photons"][str(i)]["event"])
        angle.append(readdata["photons"][str(i)]["angle"])
    stuff["pos"] = np.array(pos)
    stuff["dir"] = np.array(dir)
    stuff["wallHits"] = np.array(wallHits)
    stuff["wavelength"] = np.array(wavelength)
    stuff["event"] = np.array(event)
    stuff["angle"] = np.array(angle)

    return stuff

def angleToZ(radangles: np.ndarray):
    """radangles in radians from positive z-direction, returns angles in degrees"""
    over_pi2 = radangles[:] > np.pi/2
    # Black magic to convert 0=>180 to 0 => 90 => 0
    return radangles[:] * (180/np.pi)*(1-2*over_pi2) + 180*over_pi2

#=======================CALCULATIONS============================================================================

# Bottom of sample cell is assumed to be at z=0

def main():
    
    # HERE IS THE DATA :)
    data = readJsonData(fileName=fileName)
    
    # For selecting top or bottom exit
    botTop = sampCellZ
    if top:
        botTop = 0.

    # For printing stuff
    info = "\n"

    # Total photons
    totalPhotons = len(data["pos"])

    # pos and dir xyz
    posxyz = np.empty((totalPhotons, 3))
    dirxyz = np.empty((totalPhotons, 3))
    for i in range(totalPhotons):
        posxyz[i][0] = data["pos"][i][0]
        posxyz[i][1] = data["pos"][i][1]
        posxyz[i][2] = data["pos"][i][2]
        dirxyz[i][0] = data["dir"][i][0]
        dirxyz[i][1] = data["dir"][i][1]
        dirxyz[i][2] = data["dir"][i][2]

    # r and angle to z
    x = np.take(posxyz, 0, axis=1)
    y = np.take(posxyz, 1, axis=1)
    r = np.empty((totalPhotons))
    for i in range(len(x)):
        rr = np.sqrt(x[i]**2 + y[i]**2)
        r[i] = rr

    # Save data of photons that reached bottom or top of sample cell
    # also z position of absorbed photons
    exitPos = []
    exitR = []
    exitDir = []
    wavelenExit = []
    absZ = []
    wavelenAbs = []
    for i in range(totalPhotons):
        if data["event"][i].strip() == "exit":
            if posxyz[i][2] == botTop:
                continue
            else:
                exitPos.append(posxyz[i])
                exitDir.append(dirxyz[i])
                wavelenExit.append(data["wavelength"][i])
                exitR.append(np.sqrt(posxyz[i][0]**2 + posxyz[i][1]**2))
        elif data["event"][i].strip() == "absorption":
            absZ.append(data["pos"][i][2])
            wavelenAbs.append(data["wavelength"][i])
    exitPos = np.array(exitPos)
    exitX = np.take(exitPos, 0, axis=1)
    exitY = np.take(exitPos, 1, axis=1)
    wavelenExit = np.array(wavelenExit)
    exitDir = np.array(exitDir)
    exitR = np.array(exitR)
    absZ = np.array(absZ)
    wavelenAbs = np.array(wavelenAbs)

    # Percentage reached exit
    reachedExit = len(exitPos)
    percentage = 100 * reachedExit / totalPhotons
    info += f"{percentage:.2f} % reached exit"

    # uv and blue
    uv = 0
    blue = 0
    for i in range(len(exitPos)):
        if wavelenExit[i].strip() == "450E-9":
            blue += 1
        else:
            uv += 1
    info += f"\nUV: {uv}\nBLUE: {blue}\nCONVERTED: {100*blue/len(exitPos):.2f} %"

    # Normalize direction vectors
    normDir = np.empty((len(exitDir), 3))
    for i in range(len(exitDir)):
        length = np.sqrt(exitDir[i][0]**2 + exitDir[i][1]**2 + exitDir[i][2]**2)
        normDir[i] = [(exitDir[i][0]/length), (exitDir[i][1]/length), (exitDir[i][2]/length)]

    # Angles to plus/minus z
    ang = angleToZ(normDir)

    percentAbs = len(absZ)/totalPhotons*100
    info += f"\n{percentAbs:.2f} % of total photons were absorbed"
    absZ = np.array(absZ)
    absZax = np.linspace(0, sampCellZ, 100)
    absZBin = np.digitize(absZ, bins=absZax)
    absZperBin = np.zeros(len(absZax)+1)
    if absZBin.any():
        for point in np.nditer(absZBin):
            absZperBin[point] += 1

    #=======================PLOTS===================================================================================0

    if posDistributionPlot:
        fig1 = plt.figure(1)
        plt.hist(exitR/sampCellRadius, range=(0, 1), bins=30, color='dimgrey', rwidth=0.75)
        plt.xlabel("$r$ / R")
        plt.ylabel("Number of photons")
        if top:
            plt.title("Photons exiting sample cell at the top")
        else:
            plt.title("Photons exiting sample cell at the bottom")
        #plt.text(0.05, 19.2, f"- {percentage:.1f} % of total photons \n  reached $z$=10\n\n- of which {100*blue/len(exitPos):.1f} % blue")
        plt.tight_layout()
        if saveFigure:
            plt.savefig(posDistImName)
        plt.show()

    if angleDistributionPlot:
        fig2 = plt.figure(2)
        plt.scatter(exitR/sampCellRadius, ang, s=20, marker='.', c="mediumorchid")
        plt.xlabel("$r$ / R")
        plt.ylabel("Angle / deg")
        if top:
            plt.title("Angles of photons at sample cell top exit")
        else:
            plt.title("Angles of photons at sample cell bottom exit")
        plt.tight_layout()
        if saveFigure:
            plt.savefig(angDistImName)
        plt.show()

    if xyPlanePlot:
        # To be continued
        ymesh, xmesh = np.meshgrid(np.linspace(-1, 1, 100), np.linspace(-1, 1, 100))
        fig3, axes = plt.subplots()
        axes.scatter(exitX/sampCellRadius, exitY/sampCellRadius, marker=".", color="mediumorchid")
        axes.set_xlabel("$x$ / R")
        axes.set_ylabel("$y$ / R")
        if top:
            plt.title("Top")
        else:
            plt.title("Bottom")
        circle = Circle(xy=(0, 0), radius=1, ec="black", figure=fig3, fill=False, ls="-", visible=True, lw=2, label="$r$=R")
        axes.set_aspect(1)
        axes.add_artist(circle)
        axes.set_xlim(-2, 2)
        axes.set_ylim(-2, 2)
        plt.legend()
        #axes.pcolormesh(xmesh, ymesh, (x, y))
        plt.tight_layout()
        if saveFigure:
            plt.savefig(xyPlaneImName)
        plt.show()

    if wallHeatMapPlot:
        fig4, axes = plt.subplots()
        hm = axes.imshow(absZperBin[:,np.newaxis], cmap='jet', aspect=0.03, origin='lower')
        fig4.colorbar(hm)
        axes.set_xticks([])
        axes.set_yticks([-0.5, 100.5], ['Bot', 'Top'])
        fig4.suptitle("Photons absorbed")
        fig4.set_figwidth(2.6)
        #axes.yaxis.set_minor_locator(plt.NullLocator())
        plt.tight_layout()
        if saveFigure:
            plt.savefig(wallHeatMapImName)
        plt.show()

    if printInfo:
        print(info)

#========================================================================================================

if __name__ == "__main__":
    main()