"""Module for finding first hit object in 3D space"""

import numpy as np

import materials as mat

class Cylinder:
    """Class for a cylinder object"""
    def __init__(self, x0: np.ndarray = np.array([0,0,0]), axis: np.ndarray = np.array([0,0,0]), length: float = 1.0, radius: float = 0.2, material = mat.default, caps: np.ndarray = None):
        """
        Initialize a cylinder object
            x0: np.ndarray: Base point of the cylinder (x,y,z)
            axis: np.ndarray: Axis direction of the cylinder (3D vector)
            length: float: Length of the cylinder
            radius: float: Radius of the cylinder
            material: material: Material of the cylinder (not implemented)
            caps: np.ndarray: Boolean array for the caps of the cylinder (minrbot, maxrbot, minrtop, maxrtop)
        """
        self.x0 = x0
        self.axis = axis
        self.length = length
        self.radius = radius
        self.material = material
        if caps is None:
            self.caps = np.array([0,self.radius,0,self.radius])
        else:
            self.caps = caps
        
    def hitscan(self, start: np.ndarray, dir: np.ndarray):
        """
        Find the first intersection of a line with the cylinder
            start: np.ndarray: Starting point of the line
            dir: np.ndarray: Direction of the line
        Returns:
            (t, np.ndarray):  t value of the intersection, intersection point
        """
        # Calculate the intersection of the line with the infinite cylinder

        