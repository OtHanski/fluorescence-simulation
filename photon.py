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
                 position: np.ndarray = np.array([0,0,1]), 
                 direction: np.ndarray = None, 
                 wavelength: str = "121.567E-9"):
        # Initialize photon parameters
        self.sampCell = sampCell
        self.pos = position
        #Unit vector in random direction if no input direction
        if direction:
            self.direction = direction
        else:
            pm = np.array([-1, 1])
            direc = np.random.rand(3)
            length = np.sqrt(direc[0]**2 + direc[1]**2 + direc[2]**2)
            x1 = np.random.choice(pm) * direc[0] / length
            x2 = np.random.choice(pm) * direc[1] / length
            x3 = np.random.choice(pm) * direc[2] / length
            self.direction: np.ndarray = np.array([x1, x2, x3])
            
        self.wavelength = wavelength
        self.absorbed = False
        self.detected = False
        self.reflected = False
        self.bounces = 0

    def simulate(self):
        """Simulates photon in sampCell.
        Returns (position, direction, num of wall interactions, wavelength, exit/absorp)."""
        while True:
            hit = self.sampCell.hit_wall(position=self.pos, direction=self.direction, wavelength=self.wavelength)
            if hit[3] == "exit":
                #print("EXIT")
                return hit[0], hit[1], self.bounces, self.wavelength, hit[3]
            elif hit[3] == "specular":
                #print("SPEC")
                self.pos = hit[0]
                self.direction = hit[1]
                self.bounces += 1
            elif hit[3] == "diffuse":
                #print("DIFF")
                self.pos = hit[0]
                self.direction = hit[1]
                self.bounces += 1
            elif hit[3] == "conversion":
                #print("CONVERSION")
                self.pos = hit[0]
                self.direction = hit[1]
                self.wavelength = "450E-9"
                self.bounces += 1
            elif hit[3] == "absorption":
                #print("ABSORPTION")
                self.bounces += 1
                return hit[0], hit[1], self.bounces, self.wavelength, hit[3]
    
    def specular_reflection(self):
        pass

    def diffuse_reflection(self):
        pass

    def absorption(self):
        pass


def photonSave(fileName: str, data: list, sampCell: SampleCell = None):
    """Saves return values of simulate method (photon class) into a csv file with ';' as delimeter."""
    with open(fileName, "a") as file:
        if sampCell:
            pass #Save sample cell specs... to be implemented later
        file.write("position; direction; number of wall hits; wavelength; event")
        for point in data:
            file.write(f"\n{point[0][0], point[0][1], point[0][2]}; {point[1][0], point[1][1], point[1][2]}; {point[2]}; {point[3]}; {point[4]}")