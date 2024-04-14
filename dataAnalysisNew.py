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
plot_exitHistogram = 0
plot_angleDistribution = 1

posDistributionPlot = 1
angleDistributionPlot = 0
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
numOfWallHitsImName = "data/20240411/wallHitsHist20240411_1.png"

#===========================NOTES==========================================================================

# uv: specRefl 25%, absorp 25%, conversion 50%
# blue: specRefl 98%, absorp 2%
# 5 mm cell radius, 100 mm cell height
# 1 mm gas radius, 10 mm gas height

#========================FUNCTIONS======================================================================0===



# Read data from JSON
def readJSONData(fileName):
    """
    Reads data from JSON file and returns a dictionary with the following keys
    "cellSpecs": str,
    "pos": np.ndarray,
    "dir": np.ndarray,
    "wallHits": np.ndarray,
    "wavelength": np.ndarray,
    "event": np.ndarray,
    "angle": np.ndarray,
    "metadata": dict
    """
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

def rawJSONData(fileName):
    if not fileName:
        fileName = fh.ChooseSingleFile(initdir = "./data")
    readdata = fh.ReadJson(fileName)
    return readdata


def angleToZ(radangles: np.ndarray):
    """radangles in radians from positive z-direction, returns angles in degrees"""
    over_pi2 = radangles[:] > np.pi/2
    # Black magic to convert 0=>180 to 0 => 90 => 0
    return radangles[:] * (180/np.pi)*(1-2*over_pi2) + 180*over_pi2

def getExitRadius(data, exit = "top", wavelengths = ["450E-9", "121.567E-9"]):
    # Get photons that exited at the top, separated by wavelength
    Exits = {}
    # Set exit position
    exitZ = (exit=="top") * data["metadata"]["l_cell"]

    # If all wavelengths are wanted, return all exit radii in single array
    if wavelengths == "all":
        Exits["all"] = {}
        for key in data["photons"]:
            if data["photons"][key]["event"] == "exit" and data["photons"][key]["position"][2] == exitZ:
                Exits["all"][key] = data["photons"][key]
        Exits["all"] = np.array([Exits["all"][key]["position"][0]**2 + Exits["all"][key]["position"][1]**2 for key in Exits["all"]])
        Exits["all"] = np.sqrt(Exits["all"])
        return Exits

    for key in wavelengths:
        Exits[key] = {}
    
    for key in data["photons"]:
        if data["photons"][key]["event"] == "exit" and data["photons"][key]["position"][2] == exitZ:
            if data["photons"][key]["wavelength"] in wavelengths:
                Exits[data["photons"][key]["wavelength"]][key] = data["photons"][key]

    # Calculate the radii of the exit positions
    for wavelength in Exits:
        Exits[wavelength] = np.array([Exits[wavelength][key]["position"][0]**2 + Exits[wavelength][key]["position"][1]**2 for key in Exits[wavelength]])
        Exits[wavelength] = np.sqrt(Exits[wavelength])
    
    return Exits

def getAngles(data, exit = "top", wavelengths = "all"):
    """Returns the angles of the photons that exited the sample cell at the top or bottom, separated by wavelength
    {wavelength: [radius, angle]} or {all: [radius, angle]} if all wavelengths are wanted ("all" not finished)"""
    # Get photons that exited at the top, separated by wavelength
    Exits = {}
    # Set exit position
    exitZ = (exit=="top") * data["metadata"]["l_cell"]

    # If all wavelengths are wanted, return all exit radii in single array
    if wavelengths == "all":
        Exits["all"] = {}
        for key in data["photons"]:
            if data["photons"][key]["event"] == "exit" and data["photons"][key]["position"][2] == exitZ:
                Exits["all"][key] = data["photons"][key]
        Exits["all"] = np.array([Exits["all"][key]["position"][0]**2 + Exits["all"][key]["position"][1]**2 for key in Exits["all"]])
        Exits["all"] = np.sqrt(Exits["all"])
        return Exits

    for key in wavelengths:
        Exits[key] = {}
    
    for key in data["photons"]:
        if data["photons"][key]["event"] == "exit" and data["photons"][key]["position"][2] == exitZ:
            if data["photons"][key]["wavelength"] in wavelengths:
                Exits[data["photons"][key]["wavelength"]][key] = data["photons"][key]

    # Calculate the radii of the exit positions
    for wavelength in Exits:
        Exits[wavelength] = np.array([[Exits[wavelength][key]["position"][0]**2 + Exits[wavelength][key]["position"][1]**2, Exits[wavelength][key]["angle"]] for key in Exits[wavelength]])
        Exits[wavelength][0] = np.sqrt(Exits[wavelength][0])
        Exits[wavelength][1] = angleToZ(Exits[wavelength][1])
    
    return Exits
    
#========================PLOT FUNCTIONS=========================================================================

def posDistributionPlot(data, exit = "", blue = "450E-9", uv = "121.567E-9", savefigure = 0, logscale = False, filename = "./data/DistrHistogram.png"):
    """Plots the distribution of photons that exited the sample cell at the top or bottom as a histogram of their exit radii"""
    # Init plot
    fig = plt.figure(1)
    # Set exit position
    CellR = data["metadata"]["r_cell"]

    ExitR = getExitRadius(data, exit = exit, wavelengths = [blue, uv])
    bluer = ExitR[blue]
    uvr = ExitR[uv]

    plt.hist(bluer/CellR, range=(0, 1), bins=40, color='cornflowerblue', rwidth=0.75, alpha=0.7, label="Blue")
    plt.hist(uvr/CellR, range=(0, 1), bins=40, color="mediumorchid", rwidth=0.75, alpha=0.6, label="UV")
    plt.xlabel("$r$ / R")
    plt.ylabel("Number of photons")
    if exit == "top":
        plt.title("Photons exiting sample cell at the top")
    else:
        plt.title("Photons exiting sample cell at the bottom")
    plt.legend()
    #plt.text(0.05, 19.2, f"- {percentage:.1f} % of total photons \n  reached $z$=10\n\n- of which {100*blue/len(exitPos):.1f} % blue")
    plt.tight_layout()
    if logscale:
        plt.yscale("log")
    if savefigure:
        plt.savefig(filename)
    plt.show()


def angleDistributionPlot(data, exit = "", blue = "450E-9", uv = "121.567E-9", savefigure = 0, logscale = False, filename = "./data/DistrHistogram.png", xlim = 1):
    fig2 = plt.figure(2)
    # Set exit position
    CellR = data["metadata"]["r_cell"]
    angdata = getAngles(data, exit = exit, wavelengths = [blue, uv])

    for key in angdata:
        if key == blue:
            plt.scatter(angdata[key][:,0]/CellR, angdata[key][:,1], s=4, marker='.', c="cornflowerblue", label = "Blue")
        elif key == uv:
            plt.scatter(angdata[key][:,0]/CellR, angdata[key][:,1], s=4, marker='.', c="mediumorchid", label = "UV")
        else:
            plt.scatter(angdata[key][:,0]/CellR, angdata[key][:,1], s=4, marker='.', label = f"{key}")
    
    plt.xlabel("$r$ / R")
    plt.ylabel("Angle / deg")
    plt.xlim(0,xlim)
    if exit == "top":
        plt.title("Angles of photons at sample cell top exit")
    else:
        plt.title("Angles of photons at sample cell bottom exit")
    plt.tight_layout()
    if saveFigure:
        plt.savefig(angDistImName)
    plt.show()

#=======================CALCULATIONS============================================================================

# Bottom of sample cell is assumed to be at z=0

def main():
    
    # Read data from file
    data = rawJSONData(fileName=fileName) 
    
    # Check the metadata of the simulation
    cellexit = {"bot": 0, "top": data["metadata"]["l_cell"]}
    totalPhotons = data["metadata"]["simulations"]

    # For printing stuff
    info = "\n"

    """
    # For selecting top or bottom exit
    botTop = sampCellZ
    if top:
        botTop = 0.

    

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
    """
    #=======================PLOTS===================================================================================0
    if plot_exitHistogram:
        posDistributionPlot(data, exit = "top")
        posDistributionPlot(data, exit = "bot")
    
    if plot_angleDistribution:
        angleDistributionPlot(data, exit = "top")
        angleDistributionPlot(data, exit = "bot")

    """
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
    """

    if printInfo:
        print(info)

#========================================================================================================

if __name__ == "__main__":
    main()