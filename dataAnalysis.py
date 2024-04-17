import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from SampleCell import SampleCell
import FileHandler as fh
import time

#========================SETUP=======================================================

fileName = "./data/simulation20240417_500k.json"
if not fileName:
    fileName = fh.ChooseSingleFile(initdir = "./data")
    print(fileName)

top = 0 # Picks photons that exited at sample cell top
photonsAbsorbedLogScale = 0

# This
plot_exitHistogram = 1
# This still kinda work in progress
plot_xyPlane = 0
# Random shit not usable yet
plot_xyColorMesh = 0
# This
plot_photonsAbsorbed = 1
# This
plot_exitAngles = 1
# Not this
plot_wallHeatMap = 0
# This kinda useless
plot_photonBounces = 0
# Not this
plot_angleDistribution = 0

printInfo = 1

saveFigure = 0
posDistImName = "data/20240417/rDistribution500kBottom.png"
angDistImName = "data/20240411/simAngTop20240411_1.png"
xyPlaneImName = "data/20240411/simXYTop20240411_1.png"
wallHeatMapImName = "data/20240411/wallHeatMap20240411_1.png"
numOfWallHitsImName = "data/20240411/wallHitsHist20240411_1.png"
exitAnglesImName = "data/20240417/angleDistribution500kBottom.png"
photonsAbsorbedImName = "data/20240417/photonsAbsorbed500k.png"

#===========================NOTES==========================================================================

# uv: specRefl 25%, absorp 25%, conversion 50%
# blue: specRefl 98%, absorp 2%
# 5 mm cell radius, 100 mm cell height
# 1 mm gas radius, 10 mm gas height

#========================FUNCTIONS======================================================================0===

def saveData(fileName: str, data: list, sampCell: SampleCell = None):
    """Saves return values of simulate method (photon class) into a csv file with ';' as delimeter."""
    with open(fileName, "a") as file:
        if sampCell:
            pass # Save sample cell specs... to be implemented later
        file.write("position; direction; number of wall hits; wavelenExitgth; event")
        for point in data:
            file.write(f"\n{point[0][0], point[0][1], point[0][2]};{point[1][0], point[1][1], point[1][2]};{point[2]};{point[3]};{point[4]}")

#dat
def readDatData(fileName):
    stuff = {"cellSpecs": str, 
             "pos": np.ndarray, 
             "dir": np.ndarray, 
             "wallHits": np.ndarray, 
             "wavelength": np.ndarray,
             "event": np.ndarray}
    temp = [row.strip().split('\t') for row in open(fileName).readlines()]
    temp.remove(temp[0])
    # Clusterfuck to format data
    positions = [point[0] for point in temp]
    for i in range(len(positions)):
        positions[i] = positions[i].replace("[", "")
        positions[i] = positions[i].replace("]", "")
        positions[i] = positions[i].strip().split()
        positions[i] = [float(positions[i][0]), float(positions[i][1]), float(positions[i][2])]
    directions = [point[1] for point in temp]
    for i in range(len(directions)):
        directions[i] = directions[i].replace("[", "")
        directions[i] = directions[i].replace("]", "")
        directions[i] = directions[i].strip().split()
        directions[i] = [float(directions[i][0]), float(directions[i][1]), float(directions[i][2])]
    wallHits = [point[2] for point in temp]
    wavelengths = [point[3] for point in temp]
    events = [point[4] for point in temp]
    stuff["pos"] = np.array(positions)
    stuff["dir"] = np.array(directions)
    stuff["wallHits"] = np.array(wallHits)
    stuff["wavelength"] = np.array(wavelengths)
    stuff["event"] = np.array(events)
    
    return stuff

# Read data from JSON
def readJsonData(fileName):
    stuff = {"cellSpecs": str, 
             "pos": np.ndarray, 
             "dir": np.ndarray, 
             "wallHits": np.ndarray, 
             "wavelength": np.ndarray,
             "event": np.ndarray,
             "angle": np.ndarray,
             "metadata": dict}
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
    stuff["metadata"] = readdata["metadata"]

    return stuff

def angleToZOld(data: np.ndarray):
    if top:
        zax = np.array([0, 0, 1])
    else:
        zax = np.array([0, 0, -1])
    losAngles = []
    for a in data:
        angel = np.arccos((zax[2] * a[2]))
        losAngles.append(angel)
    return np.array(losAngles)

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
    
    # Sample cell specs -- add more if needed
    sampCellRadius = data["metadata"]["r_cell"]
    sampCellZ = data["metadata"]["l_cell"]
    

    # For selecting top or bottom exit
    botTop = sampCellZ
    if top:
        botTop = 0.

    # For printing stuff
    info = f"\n-Specs-\nCell length: {sampCellZ} m\nCell radius: {sampCellRadius} m\nWall sections: {data['metadata']["wall_sections"]}\nNumber of simulated photons: {data['metadata']["simulations"]}\nSimulation took {data['metadata']["time"]} s\n"

    # Total photons
    totalPhotons = np.size(data["event"])

    # Number of wall hits
    wallHits = np.empty((totalPhotons))
    for i in range(totalPhotons):
        wallHits[i] = int(data["wallHits"][i])
    absHits = []
    exitHits = []
    for i in range(totalPhotons):
        if data["event"][i].strip() == "absorption":
            absHits.append(int(data["wallHits"][i]))
        elif data["event"][i].strip() == "exit":
            exitHits.append(int(data["wallHits"][i]))
    absHits = np.array(absHits)
    exitHits = np.array(exitHits)


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

    # r
    x = np.take(posxyz, 0, axis=1)
    y = np.take(posxyz, 1, axis=1)
    r = np.empty((totalPhotons))
    for i in range(len(x)):
        rr = np.sqrt(x[i]**2 + y[i]**2)
        r[i] = rr
    
    # Angle to z
    angles = angleToZ(data["angle"])

    # Exit pos, r, dir, ang, wavelen of photons that reached bottom or top of sample cell
    # Separate arrays for different wavelengths
    # also z position of absorbed photons
    exitPos = []
    exitR = []
    exitDir = []
    exitAngle = []
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
                exitAngle.append(angles[i])
        elif data["event"][i].strip() == "absorption":
            absZ.append(data["pos"][i][2])
            wavelenAbs.append(data["wavelength"][i])
    # Here...
    exitPos = np.array(exitPos)
    exitX = np.take(exitPos, 0, axis=1)
    exitY = np.take(exitPos, 1, axis=1)
    exitAngle = np.array(exitAngle)

    exitXUV = []
    exitYUV = []
    exitXBlue = []
    exitYBlue = []
    exitAngleUV = []
    exitAngleBlue = []
    for i in range(len(wavelenExit)):
        if wavelenExit[i].strip() == "450E-9":
            exitXBlue.append(exitX[i])
            exitYBlue.append(exitY[i])
            exitAngleBlue.append(exitAngle[i])
        else:
            exitXUV.append(exitX[i])
            exitYUV.append(exitY[i])
            exitAngleUV.append(exitAngle[i])
    # ... And here
    wavelenExit = np.array(wavelenExit)
    exitDir = np.array(exitDir)
    exitR = np.array(exitR)
    absZ = np.array(absZ)
    wavelenAbs = np.array(wavelenAbs)
    exitXBlue = np.array(exitXBlue)
    exitYBlue = np.array(exitYBlue)
    exitXUV = np.array(exitXUV)
    exitYUV = np.array(exitYUV)
    exitAngleBlue = np.array(exitAngleBlue)
    exitAngleUV = np.array(exitAngleUV)

    absZBlue = []
    absZUV = []
    for i in range(len(wavelenAbs)):
        if wavelenAbs[i].strip() == "450E-9":
            absZBlue.append(absZ[i])
        else:
            absZUV.append(absZ[i])
    absZBlue = np.array(absZBlue)
    absZUV = np.array(absZUV)

    # Percentage reached exit...
    reachedExit = len(exitPos)
    percentage = 100 * reachedExit / totalPhotons
    if top:
        info += f"\n{percentage:.2f} % reached top exit"
    else:
        info += f"\n{percentage:.2f} % reached top exit"

    # ... of which uv and blue
    uv = []
    blue = []
    for i in range(len(exitR)):
        if wavelenExit[i].strip() == "450E-9":
            blue.append(exitR[i])
        else:
            uv.append(exitR[i])
    info += f"\nUV: {len(uv)}\nBLUE: {len(blue)}\nCONVERTED: {100*len(blue)/len(exitPos):.2f} %\n"
    uv = np.array(uv)
    blue = np.array(blue)

    # Normalize direction vectors
    normDir = np.empty((len(exitDir), 3))
    for i in range(len(exitDir)):
        length = np.sqrt(exitDir[i][0]**2 + exitDir[i][1]**2 + exitDir[i][2]**2)
        normDir[i] = [(exitDir[i][0]/length), (exitDir[i][1]/length), (exitDir[i][2]/length)]

    # Angles to plus/minus z
    #ang = angleToZ(normDir) * (180/np.pi)
    ang = angleToZ(data["angle"])

    percentAbs = len(absZ)/totalPhotons*100
    info += f"\n{percentAbs:.2f} % of total photons were absorbed\n"
    absZ = np.array(absZ)
    absZax = np.linspace(0, sampCellZ, 100)
    absZBin = np.digitize(absZ, bins=absZax)
    absZperBin = np.zeros(len(absZax)+1)
    for point in np.nditer(absZBin):
        absZperBin[point] += 1

    angleDigBins = 90
    angleDigBlue = np.digitize(exitAngleBlue, np.linspace(0, 90, angleDigBins))
    numberOfAngleBlue = np.empty(angleDigBins)
    for i in range(angleDigBins):
        numberOfAngleBlue[i] = len(np.where(angleDigBlue == i)[0])

    angleDigUV = np.digitize(exitAngleUV, np.linspace(0, 90, angleDigBins))
    numberOfAngleUV = np.empty(angleDigBins)
    for i in range(angleDigBins):
        numberOfAngleUV[i] = len(np.where(angleDigUV == i)[0])
    
    # Work in progress
    xyBins = np.empty((200, 200))

    #=======================PLOTS===================================================================================0

    # This
    if plot_exitHistogram:
        fig1 = plt.figure(1)
        plt.hist(blue/sampCellRadius, range=(0, 1), bins=40, color='skyblue', rwidth=0.8, alpha=1, label="Blue")
        plt.hist(uv/sampCellRadius, range=(0, 1), bins=40, color="purple", rwidth=0.8, alpha=0.5, label="UV")
        plt.xlabel("$r$ / R")
        plt.ylabel("Number of photons")
        if top:
            plt.title("Photon distribution at sample cell top")
        else:
            plt.title("Photon distribution at sample cell bottom")
        plt.legend()
        #plt.text(0.05, 19.2, f"- {percentage:.1f} % of total photons \n  reached $z$=10\n\n- of which {100*blue/len(exitPos):.1f} % blue")
        plt.tight_layout()
        if saveFigure:
            plt.savefig(posDistImName)
        plt.show()

    # Not this
    if plot_angleDistribution:
        fig2 = plt.figure(2)
        plt.scatter(blue/sampCellRadius, exitAngleBlue, s=4, marker='.', c="skyblue")
        plt.scatter(uv/sampCellRadius, exitAngleUV, s=4, marker='.', c="mediumorchid")
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

    # This still kinda work in progress
    if plot_xyPlane:
        fig3, (ax1, ax2) = plt.subplots(1, 2, layout="constrained", num=3)
        ax2.scatter(exitXBlue/sampCellRadius, exitYBlue/sampCellRadius, marker=".", s=1, color="skyblue")
        ax1.scatter(exitXUV/sampCellRadius, exitYUV/sampCellRadius, marker=".", s=1, color="mediumorchid")
        ax1.set_xlabel("$x$ / R")
        ax1.set_ylabel("$y$ / R")
        ax2.set_xlabel("$x$ / R")
        ax2.set_ylabel("$y$ / R")
        ax1.set_title("UV")
        ax2.set_title("Blue")
        if top:
            fig3.suptitle("Top")
        else:
            fig3.suptitle("Bottom")
        circle1 = Circle(xy=(0, 0), radius=1, ec="black", figure=fig3, fill=False, ls="-", visible=True, lw=2, label="$r$=R")
        circle2 = Circle(xy=(0, 0), radius=1, ec="black", figure=fig3, fill=False, ls="-", visible=True, lw=2, label="$r$=R")
        ax1.set_aspect(1)
        ax2.set_aspect(1)
        ax1.add_artist(circle1)
        ax2.add_artist(circle2)
        ax1.set_xlim(-1.5, 1.5)
        ax2.set_xlim(-1.5, 1.5)
        ax1.set_ylim(-1.5, 1.5)
        ax2.set_ylim(-1.5, 1.5)
        #ax1.legend()
        #ax2.legend()
        #axes.pcolormesh(xmesh, ymesh, (x, y))
        plt.tight_layout()
        if saveFigure:
            plt.savefig(xyPlaneImName)
        plt.show()

    # Random shit I wanna test
    if plot_xyColorMesh:
        fig8, ax = plt.subplots(num=8)
        z = [[exitXBlue[i], exitYBlue[i]] for i in range(len(exitYBlue))]
        ymesh, xmesh = np.meshgrid(np.linspace(-1, 1, len(exitXBlue)+1), np.linspace(-1, 1, len(exitYBlue)+1))
        ax.imshow(z, cmap="inferno", aspect="equal")
        plt.show()

    # Not this
    if plot_wallHeatMap:
        fig4, axes = plt.subplots(num=4)
        hm = axes.imshow(absZperBin[:,np.newaxis], cmap='jet', aspect=0.03, origin='lower')
        fig4.colorbar(hm)
        axes.set_xticks([])
        axes.set_yticks([-0.5, 100.5], ['Bottom', 'Top'])
        fig4.suptitle("Photons absorbed")
        fig4.set_figwidth(2.6)
        plt.tight_layout()
        if saveFigure:
            plt.savefig(wallHeatMapImName)
        plt.show()

    # This kinda useless
    if plot_photonBounces:
        # To be continued
        fig5 = plt.figure(5)
        plt.hist(absHits, bins=185, range=(0, 375), color='lightsteelblue', rwidth=1, alpha=1, label="Absorbed")
        plt.hist(exitHits, bins=185, range=(0, 375), color='darkslategrey', rwidth=1, alpha=0.5, label="Exited")
        plt.yscale("log")
        plt.xlabel("Number of wall hits")
        plt.ylabel("Number of photons")
        plt.legend()
        plt.tight_layout()
        if saveFigure:
            plt.savefig(numOfWallHitsImName)
        plt.show()

    # This
    if plot_exitAngles:
        fig7 = plt.figure(7)
        plt.plot(np.linspace(0, 90, angleDigBins), numberOfAngleBlue, color="skyblue", label="Blue")
        plt.plot(np.linspace(0, 90, angleDigBins), numberOfAngleUV, color="mediumorchid", label="UV")
        plt.xlabel("Angle / deg")
        plt.ylabel("Number of photons")
        if top:
            plt.title("Top")
        else:
            plt.title("Bottom")
        plt.legend()
        plt.tight_layout()
        if saveFigure:
            plt.savefig(exitAnglesImName)
        plt.show()

    # This
    if plot_photonsAbsorbed:
        fig4 = plt.figure(num=4)
        plt.hist(absZBlue/sampCellZ, bins=50, color='skyblue', rwidth=0.8, alpha=1, label="Blue")
        plt.hist(absZUV/sampCellZ, bins=50, color="purple", rwidth=0.8, alpha=0.5, label="UV")
        plt.title("Photons absorbed")
        plt.xlabel("$z$ / h$_{\\mathrm{cell}}$")
        plt.ylabel("Number of photons")
        plt.legend()
        plt.tight_layout()
        if photonsAbsorbedLogScale:
            plt.yscale("log")
        if saveFigure:
            plt.savefig(photonsAbsorbedImName)
        plt.show()

    if printInfo:
        print(info)

#========================================================================================================

if __name__ == "__main__":
    before = time.time()
    main()
    after = time.time()
    print(f"This shit took {after-before} seconds")