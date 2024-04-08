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

class photon:
    def __init__(self, 
                 sampCell: SampleCell,
                 position: np.ndarray = None, 
                 direction: np.ndarray = None, 
                 wavelength: str = "121.567E-9"):
        
        self.debug = 1

        if self.debug:
            print("Starting photon init")
        # Initialize photon parameters
        self.sampCell = sampCell

        # If no input position, then initialize at cell center
        if position is None:
            self.pos = sampCell.z[int(sampCell.samp/2)] 
        else:
            self.pos = position
        # Unit vector in random direction if no input direction
        self.direc = direction
        
        # In order to avoid zero vectors, generate new direction until valid
        while self.direc == None:
            self.direc = (np.random.rand(3) * 2 - 1)
            if not np.any(self.direc):
                self.direc /= np.sqrt(self.direc.dot(self.direc))
        
        self.wavelength = wavelength
        self.absorbed = False
        self.detected = False
        self.reflected = False
        self.bounces = 0
        # self.event should reflect latest event
        # "exit", "absorption", "specular", "diffuse", "conversion"
        self.event = None

        if self.debug:
            print("Photon initialized")
        
    def getDir(self):
        return self.direc

    def simulate(self):
        """Simulates photon in sampCell.
        Returns (position, direction, num of wall interactions, wavelength, exit/absorp)."""
        while True:
            hit = self.sampCell.hit_wall(position=self.pos, direction=self.direc, wavelength=self.wavelength)
            self.event = hit[3]
            if hit[3] == "exit":
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


    def data_to_string(self):
        return f"{self.pos}\t{self.direc}\t{self.bounces}\t{self.wavelength}\t{self.event}"

    def __repr__(self):
        return f"Photon at {self.pos} with direction {self.direc} and wavelength {self.wavelength}. Latest event {self.event}."
    
    def __str__(self):
        return self.__repr__()
