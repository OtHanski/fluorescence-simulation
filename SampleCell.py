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
    def __init__(self, 
                 z: np.ndarray, 
                 r: np.ndarray, 
                 specrefl: dict = {"121.567E-9": np.zeros(1000)+0.25, "450E-9": np.zeros(1000)+0.98},
                 diffrefl: dict = {"121.567E-9": np.zeros(1000), "450E-9": np.zeros(1000)},
                 absprob: dict = {"121.567E-9": np.zeros(1000)+0.25, "450E-9": np.zeros(1000)+0.02},
                 WLconversion: dict = {"121.567E-9": np.zeros(1000)+0.5, "450E-9": np.zeros(1000)},
                 wavelengths: np.ndarray = np.array(["121.567E-9", "450E-9"]),
                 samples: int = 1000
                 ):
        """
        z and r should be a numpy arrays of length samples, defining the sample cell
        wall profile as a radius r(z). 

        wavelengths should define the wavelengths being simulated

        specrefl, diffrefl, WLconversion and absprob should be dictionaries including wavelength 
        specific arrays of length samples-1, defining the specular reflection, diffuse reflection 
        and absorption probabilities as a function of z (samples-1 because we only have n-1 wall 
        pieces in a n-sample cell).
        """

        self.z = z
        self.dz = z[1] - z[0]
        self.r = r
        self.samp = samples
        self.wavelengths = wavelengths

        # Calculate the angle of the wall against the z-axis
        self.wall = np.zeros(self.samp-1)

        # Define the probabilities of reflection and absorption
        self.specular_probability = specrefl
        self.diffuse_probability = diffrefl
        self.absorption_probability = absprob
        self.WLconversion = WLconversion

        for i in range(self.samp-1):
            self.wall[i] = np.arctan((r[i+1]-r[i])/(z[i+1]-z[i]))
    
    def hit_wall(self, position = [0,0,0], direction = [1,1,1]):
        # Check location of wall hit, return exit status, new position and (specular) direction.
        # If the photon exits out of either end of the cell, return the exit location and direction.
        # Direction should be a unit vector of the form [dx,dy,dz] for easy calculations. 

        # Define transit vector for easy calculations of path.
        transit = direction / direction[2]
        
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
        pos_array = np.tile(position, (self.samp,1))
        if z_dir:
            pos_path[z_index+1:] = pos_array[z_index+1:] + transit * (self.z[z_index+1:].reshape(-1,1) - z_start)
            r_path[z_index:] =  np.sqrt(pos_path[z_index:,0]**2 + pos_path[z_index:,1]**2)
        else:
            pos_path[z_index:] = pos_array[z_index:] - transit * (self.z[z_index:].reshape(-1,1) - z_start)
            r_path[:z_index] =  np.sqrt(pos_path[:z_index,0]**2 + pos_path[:z_index,1]**2)

        # Calculate exact location of wall hit via linear interpolation
        # First, locate the index of the wall hit
        idhits = np.argwhere(np.diff(np.sign(r_path - self.r)) != 0)

        # If no wall hit, return the position and direction of the photon plus the exit status (true)
        if len(idhits) == 0:
            if z_dir:
                return [pos_path[-1,0], pos_path[-1,1], self.z[-1]], direction, True
            else:
                return [pos_path[0,0], pos_path[0,1], self.z[0]], direction, True
        # idhits returns all "hits", we only want first one.
        idhit = idhits[0][0]

        r1, r2 = r_path[idhit], r_path[idhit+1]
        w1, w2 = self.r[idhit], self.r[idhit+1]
        z1, z2 = self.z[idhit], self.z[idhit+1]
        zhit = z1 + ((w1 - r1) * (z2 - z1)) / (r2 - r1 - w2 + w1)

        poshit = position + transit * (zhit - z_start)
        surfacenormal = np.array([poshit[0], poshit[1], 0])
        refdir = direction - 2 * np.dot(direction, surfacenormal) * surfacenormal

        # Return the new position and direction, as well as the no-exit status
        return poshit, refdir, False

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

    def get_absorption_probability(self, wavelength,z):
        # Return the absorption probability for a given wavelength
        return self.absorption_probability[wavelength][self.get_z_index(z)]

    def get_specular_probability(self, wavelength, z):
        # Return the specular reflection probability for a given wavelength
        return self.specular_probability[wavelength][self.get_z_index(z)]

    def get_diffuse_probability(self, wavelength, z):
        # Return the diffuse reflection probability for a given wavelength
        return self.diffuse_probability[wavelength][self.get_z_index(z)]