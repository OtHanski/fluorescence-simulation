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
np.seterr(all = 'raise')
np.set_printoptions(precision=4)

class SampleCell:
    def __init__(self, 
                 z: np.ndarray, 
                 r: np.ndarray, 
                 specrefl: dict = None,
                 diffrefl: dict = None,
                 absprob: dict = None,
                 WLconversion: dict = None,
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

        if specrefl is None:
            specrefl = {"121.567E-9": np.zeros(samples-1)+0.25, "450E-9": np.zeros(samples-1)+0.98}
        if diffrefl is None:
            diffrefl = {"121.567E-9": np.zeros(samples-1), "450E-9": np.zeros(samples-1)}
        if absprob is None:
            absprob = {"121.567E-9": np.zeros(samples-1)+0.25, "450E-9": np.zeros(samples-1)+0.02}
        if WLconversion is None:
            WLconversion = {"121.567E-9": np.zeros(samples-1)+0.5, "450E-9": np.zeros(samples-1)}

        self.z = z
        self.dz = z[1] - z[0]
        self.r = r
        self.samp = samples
        self.wavelengths = wavelengths

        # Calculate the slope of the wall against z-axis, positive for increasing r
        self.wall = np.zeros(self.samp-1)
        # Check if perfect cylinder
        self.cylinder = not np.any(self.wall)
        

        # Define the probabilities of reflection and absorption
        self.specular_probability = specrefl
        self.diffuse_probability = diffrefl
        self.absorption_probability = absprob
        self.WLconversion = WLconversion

        for i in range(self.samp-1):
            self.wall[i] = (r[i+1]-r[i])/(z[i+1]-z[i])
    
    def hit_wall(self, position = [0,0,0], direction = [1,1,1], wavelength:str = "450E-9", verbose = 0):
        """
        Given a photon at position and direction, check if it hits the wall of the cell.

        On a hit, calculate the event type and new position/direction.

        Returns:
            np.array([x,y,z]): new position of the photon
            np.array([dx,dy,dz]): new direction of the photon
            bool: exit status
            str: event type
        """

        if self.cylinder:
            # Check if the photon hits the cylinder wall
            poshit, exitstatus = self.get_hit_location_cylinder(position, direction, verbose = verbose)
            z_index = self.get_z_index(poshit[2])# Fetch surface normal at hit location
            if z_index < 0: z_index = 0
            if z_index >= len(self.z)-1: z_index = len(self.z)-2
            surfacenormal = self.get_surfacenormal(poshit)

        else:
            # Check if the photon hits the general wall
            poshit, direction, exit, event = self.get_hit_location_general(position, direction, verbose = verbose)
            z_index = self.get_z_index(poshit[2])
        
        
        if exitstatus:
            return poshit, direction, True, "exit"

        # Resolve whether the photon reflects, absorbs or converts
        event = np.random.choice(["specular", "diffuse", "absorption", "conversion"],
                                 p = [self.specular_probability[wavelength][z_index-1], 
                                      self.diffuse_probability[wavelength][z_index-1], 
                                      self.absorption_probability[wavelength][z_index-1], 
                                      self.WLconversion[wavelength][z_index-1]])

        if event == "absorption":
            # Absorbed photon should not have a new direction
            return poshit, np.array([0,0,0]), False, event
        
        if event == "conversion":
            # Convert wavelength, random direction
            return poshit, self.randomreflection(surfacenormal), False, event

        if event == "diffuse":
            # Random direction
            return poshit,  self.randomreflection(surfacenormal), False, event
        
        if event == "specular":
            # Reflect direction specularly
            try:
                refdir = direction - 2 * np.dot(direction, surfacenormal) * surfacenormal / np.linalg.norm(surfacenormal)**2
            except Exception as er:
                raise er
            return poshit, refdir, False, event

    def get_z_index(self, z):
        # Return the index of first element before z in the z array
        
        idx = np.searchsorted(self.z, z, side="left")
        return idx-1

    def get_surfacenormal(self, pos):
        # Return the normal vector of a wall at a given position.
        x, y, z = pos
        wallslope = self.wall[self.get_z_index(z)]
        normal = np.array([-x, -y, -wallslope])
        normal = normal / np.linalg.norm(normal)
        return normal

    def randomvec(self):
        """Generates a random unit vector in 3D space. Uses spherical 
                coordinates to generate the vector to ensure uniform distribution.
        
        Returns:
            np.ndarray (x,y,z): Random unit vector"""
        theta = np.random.rand() * 2 * np.pi
        phi = (np.random.rand()-0.5) * np.pi
        z = np.sin(phi)
        x = np.cos(theta) * np.cos(phi)
        y = np.sin(theta) * np.cos(phi)
        return np.array([x,y,z])

    def randomreflection2(self, surfacenormal):
        """Return a random reflection out of a surface direction given a normal vector of said surface
        
        Returns:
            np.ndarray (x,y,z): Random reflection direction"""
        # theta is the angle in the plane of the reflection surface,
        # phi is the angle compared to normal vector (radial)
        theta = np.random.rand() * 2 * np.pi
        phi = np.random.rand() * np.pi/2

        phi2 = np.arccos(np.dot(surfacenormal, np.array([0,0,1])))

        # Calculate the new direction
        reflectionvector = np.array([np.cos(theta) * np.cos(phi) * np.sin(phi2), 
                                     np.sin(theta) * np.sin(phi) * np.cos(phi2), 
                                     np.cos(phi) * np.sin(phi2)])

        return reflectionvector
    
    def randomreflection(self, surfacenormal):
        """Return a random reflection out of a surface direction given a normal vector of said surface
        
        Returns:
            np.ndarray (x,y,z): Random reflection direction"""
        # theta is the angle in the plane of the reflection surface,
        # phi is the angle compared to normal vector (radial)
        reflectionvector = np.zeros(3)
        while not np.dot(surfacenormal, reflectionvector) > 0:
            reflectionvector = self.randomvec()

        return reflectionvector
        
    def get_hit_location_cylinder(self, startpos, dir, verbose = False):
        """Given a starting position and direction, find the location of the next wall hit in a perfect cylinder (r(z) = constant).
        
        Returns:
            np.ndarray (x,y,z): Location of wall hit
            bool: exit status"""
        
        # t^2 * (dx^2 + dy^2) + 2t(dx*x0 + dy*y0) + x0^2 + y0^2 = r^2
        a = dir[0]**2 + dir[1]**2
        b = 2*(dir[0]*startpos[0] + dir[1]*startpos[1])
        c = startpos[0]**2 + startpos[1]**2 - self.r[0]**2 # We can use r[0] as the radius is constant
        t = np.roots([a,b,c])
        # We want the positive solution
        t = t[t != 0]
        t = t[0]
        t = (-2*(dir[0]*startpos[0] + dir[1]*startpos[1]) + \
             np.sqrt(4*(dir[0]*startpos[0] + dir[1]*startpos[1])**2 - 4*(dir[0]**2 + dir[1]**2)*(startpos[0]**2 + startpos[1]**2 - self.r[0]**2)))\
             / (2*(dir[0]**2 + dir[1]**2))
        hit = startpos + t*dir

        # Check for exit status
        exitstatus =  hit[2] > self.z[-1] or hit[2] < self.z[0]
        if verbose: 
            print(f"t: {t:.4f}, hit: {hit}, startpos: {startpos}")
            print(hit, np.sqrt(hit[0]**2 + hit[1]**2))
            print(hit[2] > self.z[-1], hit[2] < self.z[0], exitstatus)
        tcorr = 0
        if exitstatus:
            exitZ = self.z[-1] if hit[2] > self.z[-1] else self.z[0]
            hit = startpos + (exitZ - startpos[2]) * dir/dir[2]
        else:
            # Check if the hit is outside the cylinder, if it is, likely floating point error => correct
            while np.sqrt(hit[0]**2 + hit[1]**2) > self.r[self.get_z_index(hit[2])]:
                # Correct the hit location
                if verbose: print("correcting")
                hit = startpos + (t-1E-8*tcorr)*dir
                tcorr += 1
        if verbose:
            print(f"Photon starting at {startpos} with direction {dir} {'hit wall' if not exitstatus else 'exited'} at {hit}, \nwith t {t:.4f}, tcorr {tcorr} and r {np.sqrt(hit[0]**2 + hit[1]**2):.4f}\n")
        
        return hit, exitstatus

    def get_hit_location_general(self, position, direction, verbose = False):
        """Old one, not working properly"""
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

        try: 
            if z_dir:
                pos_path[z_index+1:] = pos_array[z_index+1:] + transit * (self.z[z_index+1:].reshape(-1,1) - z_start)
                r_path[z_index:] =  np.sqrt(pos_path[z_index:,0]**2 + pos_path[z_index:,1]**2)
            else:
                pos_path[:z_index-1] = pos_array[:z_index-1] - transit * (self.z[:z_index-1].reshape(-1,1) - z_start)
                r_path[:z_index] =  np.sqrt(pos_path[:z_index,0]**2 + pos_path[:z_index,1]**2)
        except FloatingPointError as er:
            raise er

        # Calculate exact location of wall hit via linear interpolation
        # First, locate the index of the wall hit
        idhits = np.argwhere(np.diff(np.sign(r_path - self.r)) != 0)

        # If no wall hit, return the position and direction of the photon plus the exit status (true)
        if len(idhits) == 0:
            if z_dir:
                exitpos = position + transit * (self.z[-1] - z_start)
                return exitpos, direction, True, "exit"
            else:
                exitpos = position + transit * (self.z[0] - z_start)
                return exitpos, direction, True, "exit"
            
        # idhits returns all "hits", we only want first one.
        idhit = idhits[0][0]
        print(idhit,z_index)

        # If the photon hits the wall at the same z-index as the start, we need to handle it differently
        if idhit == z_index:
            pass
        else:
            # Linear interpolation to find exact hit location
            r1, r2 = r_path[idhit], r_path[idhit+1]
            w1, w2 = self.r[idhit], self.r[idhit+1]
            z1, z2 = self.z[idhit], self.z[idhit+1]
            zhit = z1 + ((w1 - r1) * (z2 - z1)) / ((r2 - r1) - (w2 - w1))
            print("z and rs ",self.z[idhit-2:idhit+2], r_path[idhit-2:idhit+2])

            poshit = position + transit * (zhit - z_start)
            surfacenormal = self.get_surfacenormal(poshit)

        if verbose:
            print(f"Photon starting at {position} with direction {direction} hit wall at {poshit} with surfacenormal {surfacenormal}\n\
                    Hit wall between z-indices {idhit}, {idhit+1} with r {r_path[idhit]}, {r_path[idhit+1]}\n\
                    This was between path points {pos_path[idhit]}, {pos_path[idhit+1]}")
        print(f"zhit: {zhit}, z1: {z1}, z2: {z2}, r1: {r1}, r2: {r2}, w1: {w1}, w2: {w2}\n")

    def get_absorption_probability(self, wavelength,z):
        # Return the absorption probability for a given wavelength
        return self.absorption_probability[wavelength][self.get_z_index(z)]

    def get_specular_probability(self, wavelength, z):
        # Return the specular reflection probability for a given wavelength
        return self.specular_probability[wavelength][self.get_z_index(z)]

    def get_diffuse_probability(self, wavelength, z):
        # Return the diffuse reflection probability for a given wavelength
        return self.diffuse_probability[wavelength][self.get_z_index(z)]