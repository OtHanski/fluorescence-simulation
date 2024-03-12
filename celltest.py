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

from SampleCell import SampleCell

# settings
test_cell = 1


def test_SampleCell():
    samples = 100
    z = np.linspace(0,10,samples)
    r = np.ones(samples)
    cell = SampleCell(z, r, samples = 100)

    print(f"Cell wall angle: {cell.wall}")
    print(cell.hit_wall(np.array([0,0,5]), np.array([0,0,1])))
    print(cell.hit_wall(np.array([0,0,5]), np.array([0,0,-1])))
    print(cell.hit_wall(np.array([0,0,5]), np.array([1,1,1])))

if __name__ == "__main__":
    if test_cell == 1:
        test_SampleCell()
