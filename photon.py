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

class Photon:
    def __init__(self, z, theta, phi, wavelength):
        self.z = z
        self.theta = theta
        self.phi = phi
        self.wavelength = wavelength
        self.absorption_probability = 0.5
        self.detection_probability = 0.5
        self.absorbed = False
        self.detected = False
        self.reflected = False
        self.bounces = 0

    def simulate(self, samplecell):
        # Samplecell should be a class with definition of cell geometry and reflectivity values.
        pass

    def specular_reflection(self, samplecell):
        pass

    def diffuse_reflection(self, samplecell):
        pass
