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
from photon import photon
from SampleCell import SampleCell
import FileHandler as fh

### Simulation settings ###
# Number of particles to simulate
simulations = 100
# Number of wall sections to divide the cell into (1 to n)
wall_sections = 100
# Cell parameters:
shape = "cylinder"
r_cell = 5E-3 # Radius of the cell [m]
l_cell = 100E-3 # Length of the cell [m]
# Gas cloud parameters:
gas_height = 10E-3 # Height of the gas cloud [m]
gas_offset = 5E-3 # Offset of the gas cloud from the cell bottom [m]
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


def main(simulations = simulations, wall_sections = wall_sections, r_cell = r_cell, l_cell = l_cell):
    # Generate a sample cell geometry with straight cylindrical walls
    z = np.linspace(0,l_cell,wall_sections)
    r = np.ones(wall_sections) * r_cell
    wavelengths = np.array(["121.567E-9", "450E-9"])

    # Create the parameter dictionaries for the SampleCell
    specrefl = {"121.567E-9": np.zeros(wall_sections-1)+0.25, "450E-9": np.zeros(wall_sections-1)+0.98}
    diffrefl = {"121.567E-9": np.zeros(wall_sections-1), "450E-9": np.zeros(wall_sections-1)}
    absprob = {"121.567E-9": np.zeros(wall_sections-1)+0.25, "450E-9": np.zeros(wall_sections-1)+0.02}
    WLconversion = {"121.567E-9": np.zeros(wall_sections-1)+0.5, "450E-9": np.zeros(wall_sections-1)}

    cell = SampleCell(z, r, specrefl=specrefl, diffrefl=diffrefl, absprob=absprob, 
                      WLconversion=WLconversion, samples=wall_sections, wavelengths=wavelengths)

    print(f"Beginning simulation of {simulations} photons")
    
    filepath = "./data/simulation.dat"
    with open(filepath, "w") as f:
        f.write("pos(x,y,z),dir(dx,dy,dz),bounces,wavelength,event\n")
    # Simulate the photons
    for i in range(simulations):
        # Generate a random point in the gas cloud
        if (i+1) % 100 == 0:
            print(f"Simulating photon {i+1}/{simulations}")
        pos = randomGasPoint()
        phot = photon(sampCell=cell, position=pos)
        try:
            result = phot.simulate()
            # Append result to file
            fh.WriteDat(filepath = filepath, 
                        string_to_write = phot.data_to_string()+"\n", 
                        writemode = "a")
        except:
            continue
    



# If called from CMD, run main
if __name__ == "__main__":
    main()