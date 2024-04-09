import numpy as np

def randomvec():
    """Generates a random unit vector in 3D space. Uses spherical 
    coordinates to generate the vector to ensure uniform distribution.

    Returns:
        np.ndarray (x,y,z): Random unit vector"""
    theta = np.random.rand() * 2 * np.pi
    phi = np.random.rand() * np.pi
    z = np.sin(phi)
    x = np.cos(theta) * np.cos(phi)
    y = np.sin(theta) * np.cos(phi)
    return np.array([x,y,z])

def test_randomvec():
    for i in range(10000):
        vec = randomvec()
        for j in vec:
            if j < -1 or j > 1:
                print(f"Error: {j}")