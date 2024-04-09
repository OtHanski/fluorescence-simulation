import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from SampleCell import SampleCell

#========================SETUP=======================================================

fileName = "./data/testData3.csv"
top = 1 # Picks photons that exited at sample cell top
posDistributionPlot = 0
angleDistributionPlot = 0
topDownPlot = 0
wallHeatMapPlot = 0
numOfWallHitHist = 0
printInfo = 1
posDistImName = "data/xDistribution11_11.png"
angDistImName = "data/angles1.png"
topDownImName = "data/topDowntest1.png"
heatMapImName = ""

#===========================NOTES==========================================================================

# uv: specRefl 25%, absorp 25%, conversion 50%
# blue: specRefl 98%, absorp 2%

#========================FUNCTIONS======================================================================0===

def saveData(fileName: str, data: list, sampCell: SampleCell = None):
    """Saves return values of simulate method (photon class) into a csv file with ';' as delimeter."""
    with open(fileName, "a") as file:
        if sampCell:
            pass # Save sample cell specs... to be implemented later
        file.write("position; direction; number of wall hits; wavelength; event")
        for point in data:
            file.write(f"\n{point[0][0], point[0][1], point[0][2]};{point[1][0], point[1][1], point[1][2]};{point[2]};{point[3]};{point[4]}")

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

# This is kinda useless at the moment
def plotPositions(cellRadius, posData):
    plt.hist(posData, range=(0, cellRadius/cellRadius), bins=20, color='steelblue', rwidth=0.5)

def angleToNegZ(data: np.ndarray):
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

def main():
    # HERE IS THE DATA :)
    data = readData(fileName=fileName) 

    botTop = 10.
    if top:
        botTop = 0.

    info = "\n"

    # Total photons
    totalPhotons = np.size(data["event"])

    # Save data of photons that reached bottom or top of sample cell
    exitPos = []
    wavelen = []
    for i in range(len(data["pos"])):
        if data["event"][i].strip() == "exit":
            if data["pos"][i][2] == botTop:
                continue
            else:
                exitPos.append(data["pos"][i])
                wavelen.append(data["wavelength"][i])
    exitPos = np.array(exitPos)
    wavelen = np.array(wavelen)

    # Percentage reached exit
    reachedBottom = len(exitPos)
    percentage = 100 * reachedBottom / totalPhotons
    info += f"{percentage:.2f} % reached exit"

    # uv and blue
    uv = 0
    blue = 0
    for i in range(len(exitPos)):
        if wavelen[i].strip() == "450E-9":
            blue += 1
        else:
            uv += 1
    info += f"\nUV: {uv}\nBLUE: {blue}\nCONVERTED: {100*blue/len(exitPos):.2f} %"

    # Cyl coordinate r
    x = np.take(exitPos, 0, axis=1)
    y = np.take(exitPos, 1, axis=1)
    r = []
    for i in range(len(x)):
        rr = np.sqrt(x[i]**2 + y[i]**2)
        r.append(rr)
    r = np.array(r)

    # Normalize direction vectors
    dirs = []
    for i in range(len(data["dir"])):
        if data["event"][i].strip() == "exit":
            if data["pos"][i][2] == botTop:
                continue
            else:
                length = np.sqrt(data["dir"][i][0]**2 + data["dir"][i][1]**2 + data["dir"][i][2]**2)
                xdir = data["dir"][i][0] / length
                ydir = data["dir"][i][1] / length
                zdir = data["dir"][i][2] / length
                dirlandaa = [xdir, ydir, zdir]
                dirs.append(dirlandaa)
    dirs = np.array(dirs)

    # Angles to plus/minus z
    ang = angleToNegZ(dirs) * (180/np.pi)

    # Delete bad ppoints... Throw away when debugging done
    for point in r:
        if abs(point) > 1:
            i = np.where(r==point)[0][0]
            x = np.delete(x, i)
            y = np.delete(y, i)
            r = np.delete(r, i)
            ang = np.delete(ang, i)

    # z position of absorbed photons
    absZ = []
    for i in range(len(data["event"])):
        if data["event"][i].strip() == "absorption":
            point = data["pos"][i][2]
            absZ.append(point)
    percentAbs = len(absZ)/totalPhotons*100
    info += f"\n{percentAbs:.2f} % of total photons were absorbed"
    absZ = np.array(absZ)
    absZax = np.linspace(0, 10, 100)
    absZBin = np.digitize(absZ, bins=absZax)
    absZperBin = np.zeros(len(absZax)+1)
    for point in np.nditer(absZBin):
        absZperBin[point] += 1

    #=======================PLOTS===================================================================================0

    if posDistributionPlot:
        fig1 = plt.figure(1)
        plt.hist(r, range=(0, 1), bins=30, color='dimgrey', rwidth=0.75)
        plt.xlabel("$r$ / R")
        plt.ylabel("Number of photons")
        if top:
            plt.title("Photons exiting sample cell at $z$=10")
        else:
            plt.title("Photons exiting sample cell at $z$=0")
        #plt.text(0.05, 19.2, f"- {percentage:.1f} % of total photons \n  reached $z$=10\n\n- of which {100*blue/len(exitPos):.1f} % blue")
        plt.tight_layout()
        plt.savefig(posDistImName)
        plt.show()

    if angleDistributionPlot:
        fig2 = plt.figure(2)
        plt.scatter(r, ang, s=20, marker='.', c="mediumorchid")
        plt.xlabel("$r$ / R")
        plt.ylabel("Angle / deg")
        plt.title("Angles of photons at sample cell exit")
        plt.tight_layout()
        plt.savefig(angDistImName)
        plt.show()

    if topDownPlot:
        # To be continued
        ymesh, xmesh = np.meshgrid(np.linspace(-1, 1, 100), np.linspace(-1, 1, 100))
        fig3, axes = plt.subplots()
        axes.scatter(x, y, marker=".", color="mediumorchid")
        axes.set_xlabel("$x$ / R")
        axes.set_ylabel("$y$ / R")
        circle = Circle(xy=(0, 0), radius=1, ec="black", figure=fig3, fill=False, ls="-", visible=True, lw=2)
        axes.set_aspect(1)
        axes.add_artist(circle)
        #axes.pcolormesh(xmesh, ymesh, (x, y))
        plt.tight_layout()
        plt.savefig(topDownImName)
        plt.show()

    if wallHeatMapPlot:
        fig4, axes = plt.subplots()
        hm = axes.imshow(absZperBin[:,np.newaxis], cmap='hot', aspect=0.05, origin='lower')
        fig4.colorbar(hm)
        axes.set_xticks([])
        axes.set_yticks([-0.5, 54.5, 100.5], ['Bottom', 'Sample', 'Top'])
        #axes.yaxis.set_minor_locator(plt.NullLocator())
        plt.tight_layout()
        plt.show()

    if printInfo:
        print(info)

#========================================================================================================

if __name__ == "__main__":
    main()