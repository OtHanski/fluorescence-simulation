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
import matplotlib.pyplot as plt

# Import the classes to be tested
from SampleCell import SampleCell
from photon import photon
import simulation

# Test flags
test_cell = 0
test_gas = 0
test_sim = 0
test_hits = 1

# Sample cell test settings
timetest = 1
photontest = 0

# Gas generation test settings
gas_samples = 10000
gas_height = 10E-3
gas_offset = 5E-3
gas_radius = 1E-3


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
            testphoton = photon()


def test_simulation():
    # Generate a sample cell with straight cylindrical walls (r = 1, z = 10)
    samples = 100
    z = np.linspace(0,10,samples)
    r = np.ones(samples)
    normalcell = SampleCell(z, r, samples = samples)

    # Generate fully specular cell
    specrefl = {"121.567E-9": np.zeros(samples-1)+1, "450E-9": np.zeros(samples-1)+1}
    diffrefl = {"121.567E-9": np.zeros(samples-1), "450E-9": np.zeros(samples-1)}
    absprob = {"121.567E-9": np.zeros(samples-1), "450E-9": np.zeros(samples-1)}
    WLconversion = {"121.567E-9": np.zeros(samples-1), "450E-9": np.zeros(samples-1)}
    speccell = SampleCell(z, r, specrefl=specrefl, diffrefl=diffrefl, absprob=absprob, 
                          WLconversion=WLconversion, samples=samples)

    # Generate N photons at [0,0,5], with direction from (0,0,1) to (0,0,-1), simulate the path with normal cell
    simulations = 25
    for i in range(simulations+1):
        theta = i/simulations * np.pi
        pos = np.array([0,0,5])
        direc = np.array([np.sin(theta), 0, np.cos(theta)])
        phot = photon(sampCell=normalcell, position=pos, direction=direc, id = i+1)
        phot.simulate(verbose = True)
    
    input("\n\nPress enter to continue to specular reflection\n\n")

    # Generate N photons at [0,0,5], with direction from (0,0,1) to (0,0,-1), simulate the path with specular cell
    simulations = 25
    for i in range(simulations+1):
        theta = i/simulations * np.pi
        pos = np.array([0,0,5])
        direc = np.array([np.sin(theta), 0, np.cos(theta)])
        phot = photon(sampCell=speccell, position=pos, direction=direc, id = i+1)
        phot.simulate(verbose = True)

def test_hit_wall():
    # Find the bug in the hit_wall function
    samples = 100
    z = np.linspace(0,10,samples)
    r = np.ones(samples)
    normalcell = SampleCell(z, r, samples = samples)

    # Generate fully specular cell
    specrefl = {"121.567E-9": np.zeros(samples-1)+1, "450E-9": np.zeros(samples-1)+1}
    diffrefl = {"121.567E-9": np.zeros(samples-1), "450E-9": np.zeros(samples-1)}
    absprob = {"121.567E-9": np.zeros(samples-1), "450E-9": np.zeros(samples-1)}
    WLconversion = {"121.567E-9": np.zeros(samples-1), "450E-9": np.zeros(samples-1)}
    speccell = SampleCell(z, r, specrefl=specrefl, diffrefl=diffrefl, absprob=absprob, 
                          WLconversion=WLconversion, samples=samples)

    simulations = 1
    for i in range(simulations):
        theta = 12/25 * np.pi
        pos = np.array([0,0,5])
        direc = np.array([np.sin(theta), 0, np.cos(theta)])
        phot = photon(sampCell=speccell, position=pos, direction=direc, id = i+1)
        phot.simulate(verbose = True)



def test_Gas_Generation():
    # Test the randomGasPoint function
    samples = gas_samples
    data = np.zeros((samples,3))
    rdata = np.zeros(samples)
    for i in range(samples):
        data[i] = simulation.randomGasPoint(gas_height, gas_offset, gas_radius)
        rdata[i] = np.sqrt(data[i][0]**2 + data[i][1]**2)
    
    # Plot histogram of the generated points as function of r
    plt.hist(rdata, bins=100)
    plt.show()

    # Plot the generated points in (x,y) plane with bounding circle
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.scatter(data[:,0], data[:,1], s = 0.01)
    circle = plt.Circle((0, 0), gas_radius, fill=False)
    ax.add_artist(circle)
    plt.show()



if __name__ == "__main__":
    if test_cell:
        test_SampleCell()
    if test_gas:
        test_Gas_Generation()
    if test_sim:
        test_simulation()
    if test_hits:
        test_hit_wall()
