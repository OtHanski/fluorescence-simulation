"""
Main file for running the simulation, only implement the main simulation loop here. 
Should be run from CMD.

Set up the simulation environment (Samplecell profile, how many sims), and run the 
simulation loop (Create new photon objects, simulate them). Append results to file 
every N simulations to avoid memory overflow and to avoid losing data in case of crash.

Finally, plot data and results with pyplot.

If you're *really* bored, multithread the operation :D
"""
import numpy as np
import traceback
import time

# Local modules
from photon import photon
from SampleCell import SampleCell
import FileHandler as fh

### Simulation settings ###
filename = "simulation"
output = "json" # "dat" or "json"
# Number of particles to simulate
simulations = 1000000
# Number of wall sections to divide the cell into (1 to n)
wall_sections = 150
# Cell parameters:
shape = "cylinder"
r_cell = 8E-3 # Radius of the cell [m]
l_cell = 200E-3 # Length of the cell [m]
# Gas cloud parameters:
gas_height = 5E-2 # Height of the gas cloud [m]
gas_offset = (l_cell-gas_height)/2 # Offset of the gas cloud from the cell bottom [m]
#gas_height = 0E-2 # override
gas_radius = 1E-3 # Radius of the gas cloud [m]
### End of simulation settings ###

def randomGasPoint(gas_height = gas_height, gas_offset = gas_offset, gas_radius = gas_radius):
    """
    Give a random point in the gas cloud. The gas cloud is modeled as 
    a cylinder according to the gas cloud parameters. In this function
    one should define a density distribution rho(r,z) which is used for 
    the generation.

    Returns:
        np.ndarray (x,y,z): Random point in the gas cloud, in SampleCell coordinates 
    """
    # For now we use uniform distribution
    z = np.random.rand() * gas_height + gas_offset
    r_rand = np.random.rand()
    r = np.sqrt(r_rand) * gas_radius # SQRT to get uniform distribution in volume 
    theta = np.random.rand() * 2 * np.pi

    return np.array([r*np.cos(theta), r*np.sin(theta), z])

def simulate_dat(sampCell = None, filename = filename, simulations = simulations, 
                 wall_sections = wall_sections, r_cell = r_cell, l_cell = l_cell):
    filepath = f"./data/{filename}.dat"
    with open(filepath, "w") as f:
        f.write("pos(x,y,z)\tdir(dx,dy,dz)\trad_to_z\tbounces\twavelength\tevent\n")
    # Simulate the photons
    starttime = time.time()
    writestr = ""
    for i in range(simulations):
        # Generate a random point in the gas cloud
        if (i+1) % min(int(simulations/10), 25000) == 0:
            print(f"\nSimulating photon {i+1}/{simulations}, time elapsed: {time.time()-starttime:.2f}s\n")
        pos = randomGasPoint()
        phot = photon(sampCell=sampCell, position=pos, id = i+1)
        try:
            result = phot.simulate()
            writestr += phot.data_to_string()+"\n"
            # Append results to file every 10 photons (save IO time)
            if (i+1) % 10 == 0:
                fh.WriteDat(filepath = filepath, 
                            string_to_write = writestr, 
                            writemode = "a")
                writestr = ""
        except Exception as e:
            print(f"Error at photon {i+1}: {e}")
            print(traceback.format_exc())
            continue
    print(f"Simulation of {simulations} photons completed in {time.time()-starttime:.2f}s\n")

def simulate_json(sampCell = None, filename = filename, simulations = simulations, 
                 wall_sections = wall_sections, r_cell = r_cell, l_cell = l_cell):
    filepath = f"./data/{filename}.json"
    
    # Simulate the photons
    starttime = time.time()
    data = {"photons": {},
            "metadata": {"date": time.strftime("%Y-%m-%d %H:%M:%S"),
                         "simulations": simulations,
                         "wall_sections": wall_sections,
                         "r_cell": r_cell,
                         "l_cell": l_cell}
            }
    for i in range(simulations):
        # Generate a random point in the gas cloud
        if (i+1) % min(int(simulations/10), 25000) == 0:
            print(f"\n{time.strftime('%H:%M',time.localtime())}\tSimulating photon {i+1}/{simulations}, time elapsed: {time.time()-starttime:.2f}s\n")
        pos = randomGasPoint()
        # , direction=np.array([0.02,0,1])
        phot = photon(sampCell=sampCell, position=pos, id = i+1)
        try:
            result = phot.simulate(verbose = False)
            data["photons"][i+1] = phot.data_to_dict()
        except Exception as e:
            print(f"Error at photon {i+1}: {e}")
            print(traceback.format_exc())
            continue
    data["metadata"]["time"] = time.time()-starttime
    fh.WriteJson(filepath, data)
    print(f"Simulation of {simulations} photons completed in {time.time()-starttime:.2f}s\n")


def main(simulations = simulations, wall_sections = wall_sections, r_cell = r_cell, l_cell = l_cell):
    # Generate a sample cell geometry with straight cylindrical walls
    z = np.linspace(0,l_cell,wall_sections)
    r = np.ones(wall_sections) * r_cell
    wavelengths = np.array(["121.567E-9", "450E-9"])

    # Create the parameter dictionaries for the SampleCell
    specrefl = {"121.567E-9": np.zeros(wall_sections-1)+0.10, "450E-9": np.zeros(wall_sections-1)+0.98}
    diffrefl = {"121.567E-9": np.zeros(wall_sections-1), "450E-9": np.zeros(wall_sections-1)}
    absprob = {"121.567E-9": np.zeros(wall_sections-1)+0.4, "450E-9": np.zeros(wall_sections-1)+0.02}
    WLconversion = {"121.567E-9": np.zeros(wall_sections-1)+0.5, "450E-9": np.zeros(wall_sections-1)}

    if custom:
        arrlen = len(specrefl["121.567E-9"])
        for i in range(arrlen):
            WLconv = 0.5*(i/arrlen)
            blueabs = 0.02
            uvabs = 0.25
            specrefl["121.567E-9"][i] = (1-i/arrlen)*(1-WLconv-uvabs)
            specrefl["450E-9"][i] = (1-i/arrlen)*(1-blueabs)
            diffrefl["121.567E-9"][i] = (i/arrlen)*(1-WLconv-uvabs)
            diffrefl["450E-9"][i] = (i/arrlen)*(1-blueabs)
            absprob["121.567E-9"][i] = uvabs
            absprob["450E-9"][i] = blueabs
            WLconversion["121.567E-9"][i] = WLconv
            WLconversion["450E-9"][i] = 0


    cell = SampleCell(z, r, specrefl=specrefl, diffrefl=diffrefl, absprob=absprob, 
                      WLconversion=WLconversion, samples=wall_sections, wavelengths=wavelengths)

    print(f"Beginning simulation of {simulations} photons")
    
    if output == "dat":
        simulate_dat(sampCell = cell, simulations = simulations)
    elif output == "json":
        simulate_json(sampCell = cell, simulations = simulations)
    



# If called from CMD, run main
if __name__ == "__main__":
    main()