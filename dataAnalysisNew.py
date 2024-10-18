import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import FileHandler as fh

#========================SETUP=======================================================

fileName = ""
if not fileName:
    fileName = fh.ChooseSingleFile(initdir = "./data")
    print(fileName)

plot_exitHistogram = 1
plot_angleDistribution = 1
plot_xyPlanePlot = 1
plot_WallHeatMap = 1

# Plot settings
lgscl = 1 # Logscale for histograms
densityplot = 0 # Change exithistogram to density plot

savefigures = 1
plt.rcParams['figure.figsize'] = (12, 12)
datafolder = "data/"
posfilename = "posHistogram"
angfilename = "angHistogram"
xyfilename = "xyPlanePlot"
heatmapfilename = "wallabsorb"
wallhitsfilename = "wallhits"

global dotsize

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

def getExitPosition(data, exit = "", wavelengths = ["450E-9", "121.567E-9"]):
    # Get photons that exited at the top, separated by wavelength
    Exits = {}
    # Set exit position
    exitZ = (exit == "top") * data["metadata"]["l_cell"]

    # If all wavelengths are wanted, return all exit radii in single array
    if wavelengths == "all":
        Exits["all"] = {}
        for key in data["photons"]:
            if data["photons"][key]["event"] == "exit" and data["photons"][key]["position"][2] == exitZ:
                Exits["all"][key] = data["photons"][key]
        Exits["all"] = np.array([Exits["all"][key]["position"] for key in Exits["all"]])
        return Exits

    for key in wavelengths:
        Exits[key] = {}
    
    for key in data["photons"]:
        if data["photons"][key]["event"] == "exit" and data["photons"][key]["position"][2] == exitZ:
            if data["photons"][key]["wavelength"] in wavelengths:
                Exits[data["photons"][key]["wavelength"]][key] = data["photons"][key]

    # Calculate the radii of the exit positions
    for wavelength in Exits:
        Exits[wavelength] = np.array([Exits[wavelength][key]["position"] for key in Exits[wavelength]])
    
    return Exits

def getExitRadius(data, exit = "top", wavelengths = ["450E-9", "121.567E-9"]):
    # Get photons that exited at the top, separated by wavelength
    Exits = {}
    # Set exit position
    exitZ = (exit == "top") * data["metadata"]["l_cell"]

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
        #print("WL: ", wavelength)
        Exits[wavelength] = np.array([[Exits[wavelength][key]["position"][0]**2 + Exits[wavelength][key]["position"][1]**2, Exits[wavelength][key]["angle"]] for key in Exits[wavelength]])
        #print(Exits[wavelength])
        if Exits[wavelength].size > 0:
            Exits[wavelength][:,0] = np.sqrt(Exits[wavelength][:,0])
            Exits[wavelength][:,1] = angleToZ(Exits[wavelength][:,1])
    
    return Exits

def getAbsorbZ(data, wavelengths = "all"):
    """Returns the z position of the photons that were absorbed in the sample cell, separated by wavelength"""
    Absorbed = {}
    if wavelengths == "all":
        Absorbed["all"] = {}
        for key in data["photons"]:
            if data["photons"][key]["event"] == "absorption":
                Absorbed["all"][key] = data["photons"][key]
        Absorbed["all"] = np.array([Absorbed["all"][key]["position"][2] for key in Absorbed["all"]])
        return Absorbed
    else:
        for wavelength in wavelengths:
            Absorbed[wavelength] = {}
        for key in data["photons"]:
            if data["photons"][key]["event"] == "absorption" and (data["photons"][key]["wavelength"] in wavelengths):
                Absorbed[data["photons"][key]["wavelength"]][key] = data["photons"][key]
        for wavelength in Absorbed:
            Absorbed[wavelength] = np.array([Absorbed[wavelength][key]["position"][2] for key in Absorbed[wavelength]])
        return Absorbed

#========================PLOT FUNCTIONS=========================================================================

def posDistributionPlot(data, exit = "", blue = "450E-9", uv = "121.567E-9", savefigure = 0, logscale = True, filename = "./data/DistrHistogram.png"):
    """Plots the distribution of photons that exited the sample cell at the top or bottom as a histogram of their exit radii"""
    # Init plot
    fig = plt.figure(1)
    # Set exit position
    CellR = data["metadata"]["r_cell"]
    # Number of bins
    binN = 40

    ExitR = getExitRadius(data, exit = exit, wavelengths = [blue, uv])
    bluer = ExitR[blue]
    uvr = ExitR[uv]
    weightsblue = np.zeros(len(bluer))+1
    weightsuv = np.zeros(len(uvr))+1
    # If calculating density, normalise by 2*pi*r
    if densityplot:
        weightsblue = 1/(2*np.pi*bluer)
        weightsuv = 1/(2*np.pi*uvr)

    plt.hist(bluer/CellR, range=(0, 1), bins=binN, weights = weightsblue, color='cornflowerblue', rwidth=0.75, alpha=0.7, label="Blue")
    plt.hist(uvr/CellR, range=(0, 1), bins=binN, weights = weightsuv, color="mediumorchid", rwidth=0.75, alpha=0.6, label="UV")
    plt.xlabel("$r$ / R")
    plt.ylabel(f"Number of photons{' per unit area' if densityplot else ''}")
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

def angleDistributionPlot(data, exit = "", blue = "450E-9", uv = "121.567E-9", savefigure = 0, filename = "./data/AngleHistogram.png", xlim = 1):
    fig2 = plt.figure(2)
    # Set exit position
    CellR = data["metadata"]["r_cell"]
    angdata = getAngles(data, exit = exit, wavelengths = [blue, uv])

    for key in angdata:
        if angdata[key].size > 0:
            if key == blue:
                plt.scatter(angdata[key][:,0]/CellR, angdata[key][:,1], s=dotsize, marker='.', c="cornflowerblue", label = "Blue")
            elif key == uv:
                plt.scatter(angdata[key][:,0]/CellR, angdata[key][:,1], s=dotsize, marker='.', c="mediumorchid", label = "UV")
            else:
                plt.scatter(angdata[key][:,0]/CellR, angdata[key][:,1], s=dotsize, marker='.', label = f"{key}")
        
    plt.xlabel("$r$ / R")
    plt.ylabel("Angle / deg")
    plt.xlim(0,xlim)
    if exit == "top":
        plt.title("Angles of photons at sample cell top exit")
    else:
        plt.title("Angles of photons at sample cell bottom exit")
    plt.tight_layout()
    if savefigure:
        plt.savefig(filename)
    plt.show()

def xyPlanePlot(data, exit = "", blue = "450E-9", uv = "121.567E-9", savefigure = 0, filename = "./data/xyPlanePlot.png"):
    fig3, axes = plt.subplots(num=3)

    # Set exit position
    CellR = data["metadata"]["r_cell"]
    exitPos = getExitPosition(data, exit = exit, wavelengths = [blue, uv])

    if exitPos[blue].size > 0:
        axes.scatter(exitPos[blue][:,0]/CellR, exitPos[blue][:,1]/CellR, marker=".", s=dotsize, color="skyblue", label = "Blue")
    if exitPos[uv].size > 0:
        axes.scatter(exitPos[uv][:,0]/CellR, exitPos[uv][:,1]/CellR, marker=".", s=dotsize, color="mediumorchid", label = "UV")
    axes.set_xlabel("$x$ / R")
    axes.set_ylabel("$y$ / R")
    if exit == "top":
        plt.title("Top")
    elif exit == "bot":
        plt.title("Bottom")
    circle = Circle(xy=(0, 0), radius=1, ec="black", figure=fig3, fill=False, ls="-", visible=True, lw=2, label="$r$=R")
    axes.set_aspect(1)
    axes.add_artist(circle)
    axes.set_xlim(-2, 2)
    axes.set_ylim(-2, 2)
    plt.legend()
    #axes.pcolormesh(xmesh, ymesh, (x, y))
    plt.tight_layout()
    if savefigure:
        plt.savefig(filename)
    plt.show()

def wallHeatMapPlot(data, bins = 100, blue = "450E-9", uv = "121.567E-9", logscale = 1, savefigure = 0, filename = "./data/HeatMapPlot.png"):
    fig4 = plt.figure(num=4)

    Z_absorption = getAbsorbZ(data, wavelengths = [blue, uv])
    plt.hist(Z_absorption[blue], bins=bins, color='cornflowerblue', rwidth=0.75, alpha=0.7, label="Blue")
    plt.hist(Z_absorption[uv], bins=bins, color="mediumorchid", rwidth=0.75, alpha=0.6, label="UV")

    plt.title("Photons absorbed")
    plt.tight_layout()
    if logscale:
        plt.yscale("log")
    if savefigure:
        plt.savefig(filename)
    plt.show()

#=======================CALCULATIONS============================================================================

# Bottom of sample cell is assumed to be at z=0

def main(fileName = fileName, plot_exitHistogram = plot_exitHistogram, plot_angleDistribution = plot_angleDistribution, 
         plot_xyPlanePlot = plot_xyPlanePlot, plot_WallHeatMap = plot_WallHeatMap):
    
    # Read data from file
    data = rawJSONData(fileName=fileName)

    # Set dotsize according to sample size
    global dotsize 
    dotsize = 16*(20000/data["metadata"]["simulations"])**0.85 # Exponent is arbitrary, this one just looks good
    print("Dotsize: ", dotsize)

    #=======================PLOTS===================================================================================0
    if plot_exitHistogram:
        posDistributionPlot(data, exit = "top", savefigure=savefigures, filename=datafolder+posfilename+"Top.png", logscale = lgscl)
        posDistributionPlot(data, exit = "bot", savefigure=savefigures, filename=datafolder+posfilename+"Bot.png", logscale = lgscl)
    
    if plot_angleDistribution:
        angleDistributionPlot(data, exit = "top", savefigure=savefigures, filename=datafolder+angfilename+"Top.png")
        angleDistributionPlot(data, exit = "bot", savefigure=savefigures, filename=datafolder+angfilename+"Bot.png")
    
    if plot_xyPlanePlot:
        xyPlanePlot(data, exit = "top", savefigure=savefigures, filename=datafolder+xyfilename+"Top.png")
        xyPlanePlot(data, exit = "bot", savefigure=savefigures, filename=datafolder+xyfilename+"Bot.png")
    
    if plot_WallHeatMap:
        wallHeatMapPlot(data, savefigure=savefigures, filename=datafolder+heatmapfilename+".png")

#========================================================================================================

if __name__ == "__main__":
    main()