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

class SampleCell:
    def __init__(self, z, r, z_samples, r_samples):
        self.z = z
        self.r = r
        self.zlen = z_samples
        self.wall = np.zeros(z_samples-1)

        for i in range(z_samples-1):
            self.wall[i] = np.arctan((r[i+1]-r[i])/(z[i+1]-z[i]))
    
    def hit_wall(self, z_start, r_start, theta, phi):
        # Check location of wall hit, return new z, r, theta, phi
        pass

        
    def get_z_index(self, z):
        # Return the index of z in the z array
        pass
        """
        def find_nearest(array,value):
            idx = np.searchsorted(array, value, side="left")
            if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
                return array[idx-1]
            else:
                return array[idx]
        """

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