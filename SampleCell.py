"""
File for defining the sample cell geometry. The cell has a 
z-coordinate and a radius as a function of z. The profile should
be defined as a np-array of length n of r values corresponding to z values. 
Constructor of the form SampleCell(samples, z, r), with z and r numpy arrays 
of length samples.

Then calculate the angle of the wall against z-axis as a function of z, defined 
as the difference between successive r values (angle = atan(dr/dz). Save these 
angles in a numpy array of length n-1 for quick lookup during photon reflection events.

For the same buckets, define a wavelength-dependent absorption probability on the wall
as well as probabilities of specular and diffuse reflection.
"""

class SampleCell:
    def __init__(self, samples, z, r):
        pass

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