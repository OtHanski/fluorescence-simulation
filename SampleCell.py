"""
File for defining the sample cell geometry. The cell has a 
z-coordinate and a radius as a function of z. We always assume
cylindrical symmetry. The profile should
be defined as a np-array of length n of r values corresponding to z values. 
Constructor of the form SampleCell(samples, z, r), with z and r numpy arrays 
of length samples.

Then calculate the angle of the wall against z-axis as a function of z, defined 
as the difference between successive r values (angle = atan(dr/dz). Save these 
angles in a numpy array of length n-1 for quick lookup during photon reflection events.

For the same buckets, define a wavelength-dependent absorption probability on the wall
as well as probabilities of specular and diffuse reflection.

To account for non-cylindrical (but still cylindrically symmetric) cells, we must check 
the wall hits in 3D, as hits to an angled wall will change the direction of the photon 
for both angles.
"""

import numpy as np
import math

class SampleCell:
    def __init__(self, z, r, samples = 1000):
        self.z = z
        self.dz = z[1] - z[0]
        self.r = r
        self.samp = samples

        # Calculate the angle of the wall against the z-axis
        self.wall = np.zeros(self.samp-1)

        for i in range(self.samp-1):
            self.wall[i] = np.arctan((r[i+1]-r[i])/(z[i+1]-z[i]))
    
    def hit_wall(self, position = [0,0,0], direction = [1,1,1]):
        # Check location of wall hit, return exit status, new position and (specular) direction.
        # If the photon exits out of either end of the cell, return the exit location and direction.
        # Direction should be of the form [dx,dy,1] for easy calculations. If not, then we correct.

        direction = direction / direction[2]
        
        z_start = position[2]
        # Check direction of photon
        if direction[2] > 0:
            z_dir = True
        else:
            z_dir = False
        # Find the location of z_start in the z array:
        z_index = self.get_z_index(z_start)

        # Construct the r coordinates of the photon path
        r_path = np.zeros((self.samp))
        pos_path = np.zeros((self.samp,3))
        if z_dir:
            pos_path[z_index:] = position + direction * (self.z[z_index:] - z_start)
            r_path[z_index:] =  np.sqrt(pos_path[z_index:,0]**2 + pos_path[z_index:,1]**2)
        else:
            pos_path[:z_index] = position - direction * (z_start - self.z[:z_index])
            r_path[:z_index] =  np.sqrt(pos_path[:z_index,0]**2 + pos_path[:z_index,1]**2)

        # Calculate exact location of wall hit via linear interpolation
        idhit = np.argwhere(np.diff(np.sign(r_path - self.r)) != 0)

        r1, r2 = r_path[idhit], r_path[idhit+1]
        r3, r4 = self.r[idhit], self.r[idhit+1]
        z1 = self.z[idhit]
        zhit = z1 + (r3 - r1) / (r2-r4 + r3-r1)

        poshit = position + direction * (zhit - z_start)
        
        

        
    def get_z_index(self, z):
        # Return the index of z in the z array
        
        idx = np.searchsorted(self.z, z, side="left")
        if idx > 0 and (idx == len(self.z) or math.fabs(z - self.z[idx-1]) < math.fabs(z - self.z[idx])):
            return idx-1
        else:
            return idx

    def get_wall_angle(self, z):
        # Return the angle of the wall at z
        pass 

    def get_absorption_probability(self, wavelength):
        # Return the absorption probability for a given wavelength
        return self.absorption_probability

    def get_specular_probability(self, wavelength):
        # Return the specular reflection probability for a given wavelength
        return self.specular_probability

    def get_diffuse_probability(self, wavelength):
        # Return the diffuse reflection probability for a given wavelength
        return self.diffuse_probability