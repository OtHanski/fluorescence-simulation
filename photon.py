"""
File for the photon class. Photon should at least have a z-coordinate, 
two angles indicating direction, and a wavelength. Define also a wavelength-dependent
absorption probability on the wall and detection probability on the SiPM.

Simulation of particle should be run as method of object. In our case we will use a 
cylindrically symmetric system => We can calculate bounces in the z-plane and then 
use the angle to calculate the number of bounces in the r-plane.

For constructor, create a random angle from cell center to the wall (with some
freely adjustable distribution function), calculate initial z-position from this.
Also create z-min/z-max to model length of sample cell.

For UV photons, ~50% chance of TPB conversion (check papers for accurate value),
otherwise handle as reflection event. When converted, change wavelength and 
randomize angles in both planes.

For reflecion event, check Samplecell for angle of wall at current z,
check for specular or diffuse reflection (and absorption if neither), and calculate 
new angles accordingly. Then move photon to new z position and repeat until either 
end of cell is reached or photon absorbed.

Upon hitting either end of the tube, record exit position and angle.
"""
from SampleCell import SampleCell
import numpy as np
import traceback

class photon:
    def __init__(self, 
                 sampCell: SampleCell,
                 position: np.ndarray = np.array([np.nan,np.nan,np.nan]), 
                 direction: np.ndarray = np.array([0,0,0]), 
                 wavelength: str = "121.567E-9",
                 id = 0):
        # Debug flags
        self.debug = 0
        self.logold = 0

        # Initialize photon parameters
        self.sampCell = sampCell
        self.id = id

        # If no input position, then initialize at cell center
        if np.nan in position:
            self.pos = sampCell.z[int(sampCell.samp/2)] 
        else:
            self.pos = position

        self.direc = direction
        # Unit vector in random direction if no input direction
        if np.all(self.direc == 0):
            # Use the sample cell function for proper random gen
            self.direc = self.sampCell.randomvec()
        
        self.wavelength = wavelength
        self.bounces = 0
        # self.event should reflect latest event
        # "exit", "absorption", "specular", "diffuse", "conversion"
        self.event = None
        
    def getDir(self):
        return self.direc

    def simulate(self):
        """Simulates photon in sampCell.
        Returns (position, direction, num of wall interactions, wavelength, exit/absorp)."""
        try:
            while True:
                if self.logold:
                    self.oldevent = self.event
                    self.oldpos = self.pos
                    self.olddirec = self.direc

                hit = self.sampCell.hit_wall(position=self.pos, direction=self.direc, 
                                            wavelength=self.wavelength, debug = self.debug)
                self.event = hit[3]
                if hit[3] == "exit":
                    """ Commented out to avoid excessive printouts, add back if problems arise
                    if self.debug: 
                        print(f"Photon exited at {hit[0]} with direction {hit[1]}, old event {self.oldevent}, old pos {self.oldpos}, old direc {self.olddirec}")
                    if ((hit[0] == np.array([0,0,0])).all() or (hit[0] == np.array([0,0,10])).all()): 
                        print(f"Photon exited at {hit[0]} with direction {hit[1]}, \nold event {self.oldevent}, old pos {self.oldpos}, old direc {self.olddirec}")
                    """
                    self.pos = hit[0]
                    self.direc = hit[1]
                    return hit[0], hit[1], self.bounces, self.wavelength, hit[3]
                elif hit[3] == "specular":
                    self.pos = hit[0]
                    self.direc = hit[1]
                    self.bounces += 1
                elif hit[3] == "diffuse":
                    self.pos = hit[0]
                    self.direc = hit[1]
                    self.bounces += 1
                elif hit[3] == "conversion":
                    self.pos = hit[0]
                    self.direc = hit[1]
                    self.wavelength = "450E-9"
                    self.bounces += 1
                elif hit[3] == "absorption":
                    self.pos = hit[0]
                    self.direc = hit[1]
                    return hit[0], hit[1], self.bounces, self.wavelength, hit[3]
        except Exception as e:
            print(traceback.format_exc())
            print(self)
            return None


    def data_to_string(self):
        return f"{self.pos}\t{self.direc}\t{self.bounces}\t{self.wavelength}\t{self.event}"

    def __repr__(self):
        return f"Photon {self.id} at {self.pos} with direction {self.direc} and wavelength {self.wavelength}. Latest event {self.event}."
    
    def __str__(self):
        return self.__repr__()
