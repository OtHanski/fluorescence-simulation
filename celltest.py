"""
Test the SampleCell class:
- Generate a sample cell with straight cylindrical walls (r = 1, z = 10)
    - Print the wall angles (should be 0)
    - Generate photons at the center of the cell, and check if it hits the wall
        - Direction should be [0,0,1], so it should exit along z-axis, positive dir
        - Direction should be [0,0,-1], so it should exit along z-axis, negative dir
        - Direction should be [1,1,1], so it should hit the wall, be reflected at [-1,-1,1]
        - Direction should be [1,2,2], so it should hit the wall, be reflected at [-1,-2,2]
- Test the execution time of the hit_wall function
- Generate a photon and test the photon class
    - Generate 100 photons at [0,0,5], with random direction, simulate the path
    - Plot the exit positions and directions
"""

import numpy as np
import random as rng
import time

# Import the classes to be tested
from SampleCell import SampleCell
from photon import Photon

# settings
test_cell = 1
timetest = 1
photontest = 0


def test_SampleCell():

    # Generate a sample cell with straight cylindrical walls (r = 1, z = 10)
    samples = 1000
    z = np.linspace(0,10,samples)
    r = np.ones(samples)
    cell = SampleCell(z, r, samples = samples)

    # Test execution time of hit_wall function with 10k iterations 
    # Seems to be O(samples) exec time, 0.45 ms for 10k samples, 5 ms for 100k samples
    if timetest:
        print("Testing execution time of hit_wall function")
        start_time = time.time()
        for i in range(10000):
            cell.hit_wall(np.array([0,0,5]), np.array([rng.random(),rng.random(),rng.random()]))
        end_time = time.time()
        print(f"Execution time: {(end_time - start_time)/10:.10f} ms")

    # Plain wall hit tests
    print("\nTesting wall hits, generated at 0,0,5")
    test1 = cell.hit_wall(np.array([0,0,5]), np.array([0,0,1]))
    print(f"\nTest for direction [0,0,1]\nHit location: {test1[0]}, \ndirection: {test1[1]}, \nexit status: {test1[2]}")
    test2 = cell.hit_wall(np.array([0,0,5]), np.array([0,0,-1]))
    print(f"\nTest for direction [0,0,-1]\nHit location: {test2[0]}, \ndirection: {test2[1]}, \nexit status: {test2[2]}")
    test3 = cell.hit_wall(np.array([0,0,5]), np.array([1,1,1]))
    print(f"\nTest for direction [1,1,1]\nHit location: {test3[0]}, \ndirection: {test3[1]}, \nexit status: {test3[2]}")
    test4 = cell.hit_wall(np.array([0,0,5]), np.array([1,2,2]))
    print(f"\nTest for direction [1,2,2]\nHit location: {test4[0]}, \ndirection: {test4[1]}, \nexit status: {test4[2]}")

    # Test the photon class
    if photontest:
        # TODO later, when photon class is implemented
        for i in range(100):
            photon = Photon()

if __name__ == "__main__":
    if test_cell == 1:
        test_SampleCell()
