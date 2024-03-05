"""
Test the SampleCell class:
- Generate a sample cell with straight cylindrical walls (r = 1, z = 10)
    - Print the wall angles (should be 0)
    - Generate photons at the center of the cell, and check if it hits the wall
        - Direction should be [0,0,1], so it should exit along z-axis, positive dir
        - Direction should be [0,0,-1], so it should exit along z-axis, negative dir
        - Direction should be [1,1,1], so it should hit the wall, be reflected at [-1,-1,1]
"""

import numpy as np
import math

from ../SampleCell import SampleCell

def test_SampleCell():
    z = np.linspace(0,10,100)
    r = np.ones(100)
    cell = SampleCell(z, r)

    print(cell.wall)
    print(cell.hit_wall([0,0,5], [0,0,1]))
    print(cell.hit_wall([0,0,5], [0,0,-1]))
    print(cell.hit_wall([0,0,5], [1,1,1]))
