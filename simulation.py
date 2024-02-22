"""
Main file for running the simulation, only implement the main simulation loop here. 
Should be run from CMD.

Set up the simulation environment (Samplecell profile, how many sims), and run the 
simulation loop (Create new photon objects, simulate them). Append results to file 
every N simulations to avoid memory overflow and to avoid losing data in case of crash.

Finally, plot data and results with pyplot.

If you're *really* bored, multithread the operation :D
"""

from photon import Photon
from SampleCell import SampleCell

### Simulation settings ###
# Number of particles to simulate
simulations = 1000

# Number of wall sections to divide the cell into (1 to n)
wall_sections = 100



def main():
    pass






# If called from CMD, run main
if __name__ == "__main__":
    main()