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
posDistImName = "data/20240411/simDistTop20240411_1.png"
angDistImName = "data/20240411/simAngTop20240411_1.png"
xyPlaneImName = "data/20240411/simXYTop20240411_1.png"
wallHeatMapImName = "data/20240411/wallHeatMap20240411_1.png"
numOfWallHitsImName = "data/20240411/wallHitsHist20240411_1.png"

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

# Old
def readData(fileName) -> dict:
    """First row in file should include sample cell specs. 
    Second row is ignored, and the rest is saved to a list as tuples."""
    stuff = {"cellSpecs": str, 
             "pos": np.ndarray, 
             "dir": np.ndarray, 
             "wallHits": np.ndarray, 
             "wavelength": np.ndarray,
             "event": np.ndarray}
    with open(fileName) as file:
        tempdata = []
        for row in file:
            rdata = row.strip().split(';')
            tempdata.append(rdata)
    tempdata.remove(tempdata[0])
    positions = []
    directions = []
    wallHits = []
    wavelengths = []
    events = []
    # Throw away when debugging done
    for point in tempdata:
        try:
            positions.append(eval(point[0]))
        except:
            continue
        try:
            directions.append(eval(point[1]))
        except:
            positions.remove(positions[-1])
            continue
        try:
            wallHits.append(int(point[2]))
        except:
            positions.remove(positions[-1])
            directions.remove(directions[-1])
            continue
        try:
            wavelengths.append(point[3])
        except:
            positions.remove(positions[-1])
            directions.remove(directions[-1])
            wallHits.remove(wallHits[-1])
            continue
        try:
            events.append(point[4])
        except:
            positions.remove(positions[-1])
            directions.remove(directions[-1])
            wallHits.remove(wallHits[-1])
            wavelengths.remove(wavelengths[-1])
            continue
    stuff["pos"] = np.array(positions)
    stuff["dir"] = np.array(directions)
    stuff["wallHits"] = np.array(wallHits)
    stuff["wavelength"] = np.array(wavelengths)
    stuff["event"] = np.array(events)
    return stuff

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

def angleToZ(data: np.ndarray):
    if top:
        zax = np.array([0, 0, 1])
    else:
        zax = np.array([0, 0, -1])
    losAngles = []
    for a in data:
        angel = np.arccos((zax[2] * a[2]))
        losAngles.append(angel)
    return np.array(losAngles)

#=======================CALCULATIONS============================================================================

# Bottom of sample cell is assumed to be at z=0

def main():
    
    # HERE IS THE DATA :)
    data = readDatData(fileName=fileName) 
    
    # For selecting top or bottom exit
    botTop = sampCellZ
    if top:
        botTop = 0.

    # For printing stuff
    info = "\n"

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
    exitXUV = []
    exitYUV = []
    exitXBlue = []
    exitYBlue = []
    for i in range(len(wavelenExit)):
        if wavelenExit[i].strip() == "450E-9":
            exitXBlue.append(exitX[i])
            exitYBlue.append(exitY[i])
        else:
            exitXUV.append(exitX[i])
            exitYUV.append(exitY[i])
    wavelenExit = np.array(wavelenExit)
    exitDir = np.array(exitDir)
    exitR = np.array(exitR)
    absZ = np.array(absZ)
    wavelenAbs = np.array(wavelenAbs)
    exitXBlue = np.array(exitXBlue)
    exitYBlue = np.array(exitYBlue)
    exitXUV = np.array(exitXUV)
    exitYUV = np.array(exitYUV)

    # Percentage reached exit
    reachedExit = len(exitPos)
    percentage = 100 * reachedExit / totalPhotons
    if top:
        info += f"{percentage:.2f} % reached top exit\n"
    else:
        info += f"{percentage:.2f} % reached top exit\n"

    # uv and blue
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
    ang = angleToZ(normDir) * (180/np.pi)

    percentAbs = len(absZ)/totalPhotons*100
    info += f"\n{percentAbs:.2f} % of total photons were absorbed\n"
    absZ = np.array(absZ)
    absZax = np.linspace(0, sampCellZ, 100)
    absZBin = np.digitize(absZ, bins=absZax)
    absZperBin = np.zeros(len(absZax)+1)
    for point in np.nditer(absZBin):
        absZperBin[point] += 1

    #=======================PLOTS===================================================================================0

    if posDistributionPlot:
        fig1 = plt.figure(1)
        plt.hist(blue/sampCellRadius, range=(0, 1), bins=40, color='cornflowerblue', rwidth=0.75, alpha=0.7, label="Blue")
        plt.hist(uv/sampCellRadius, range=(0, 1), bins=40, color="mediumorchid", rwidth=0.75, alpha=0.6, label="UV")
        plt.xlabel("$r$ / R")
        plt.ylabel("Number of photons")
        if top:
            plt.title("Photons exiting sample cell at the top")
        else:
            plt.title("Photons exiting sample cell at the bottom")
        plt.legend()
        #plt.text(0.05, 19.2, f"- {percentage:.1f} % of total photons \n  reached $z$=10\n\n- of which {100*blue/len(exitPos):.1f} % blue")
        plt.tight_layout()
        if saveFigure:
            plt.savefig(posDistImName)
        plt.show()

    if angleDistributionPlot:
        fig2 = plt.figure(2)
        plt.scatter(exitR/sampCellRadius, ang, s=4, marker='.', c="mediumorchid")
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
        fig3, axes = plt.subplots(num=3)
        axes.scatter(exitXBlue/sampCellRadius, exitYBlue/sampCellRadius, marker=".", s=4, color="skyblue")
        axes.scatter(exitXUV/sampCellRadius, exitYUV/sampCellRadius, marker=".", s=4, color="mediumorchid")
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
        fig4, axes = plt.subplots(num=4)
        hm = axes.imshow(absZperBin[:,np.newaxis], cmap='jet', aspect=0.03, origin='lower')
        fig4.colorbar(hm)
        axes.set_xticks([])
        axes.set_yticks([-0.5, 100.5], ['Bot', 'Top'])
        fig4.suptitle("Photons absorbed")
        fig4.set_figwidth(2.6)
        plt.tight_layout()
        if saveFigure:
            plt.savefig(wallHeatMapImName)
        plt.show()

    if numOfWallHitHist:
        # To be continued
        fig5 = plt.figure(5)
        plt.hist(absHits, bins=101, range=(0, 100), color='cornflowerblue', rwidth=0.75, alpha=0.5, align="mid", label="Absorbed")
        plt.hist(exitHits, bins=101, range=(0, 100), color='mediumorchid', rwidth=0.75, alpha=0.5, label="Exited")
        plt.legend()
        plt.tight_layout()
        if saveFigure:
            plt.savefig(numOfWallHitsImName)
        plt.show()

    if printInfo:
        print(info)

#========================================================================================================

if __name__ == "__main__":
    main()