import numpy as np

def random_unit_vector(rng=np.random):
    """Return a random unit vector uniformly distributed on the sphere."""
    phi = rng.uniform(0, 2 * np.pi)
    costheta = rng.uniform(-1, 1)
    sintheta = np.sqrt(1 - costheta ** 2)
    return np.array([sintheta * np.cos(phi), sintheta * np.sin(phi), costheta])

class Photon:
	def __init__(self, position, direction, wavelength, absorbed=False):
		"""
		position: (x, y, z) tuple or array
		direction: (dx, dy, dz) tuple or array (should be normalized)
		wavelength: string, e.g. 'blue' or 'uv'
		absorbed: bool, True if photon has been absorbed
		"""
		self.position = np.array(position, dtype=float)
		self.direction = np.array(direction, dtype=float) / np.linalg.norm(direction)
		self.wavelength = wavelength
		self.absorbed = absorbed

	def move(self, distance):
		"""Move the photon along its direction by the given distance."""
		self.position += self.direction * distance

	def absorb(self):
		"""Mark the photon as absorbed."""
		self.absorbed = True

	def convert_wavelength(self):
		"""Convert the photon's wavelength (e.g., from 'uv' to 'blue')."""
		if self.wavelength == 'uv':
			self.wavelength = 'blue'
		# else: no conversion

	def specular_reflect(self, surface, cylinder):
		"""Reflect the photon specularly off the surface."""
		n = None
		if surface == 'side':
			# Normal is radial direction
			rel = self.position - cylinder.center
			h = np.dot(rel, cylinder.axis)
			proj = cylinder.center + h * cylinder.axis
			n = (self.position - proj)
			n = n / np.linalg.norm(n)
		elif surface == 'top':
			n = cylinder.axis
		elif surface == 'bottom':
			n = -cylinder.axis
		else:
			raise ValueError('Unknown surface for reflection')
		self.direction = self.direction - 2 * np.dot(self.direction, n) * n
		self.direction /= np.linalg.norm(self.direction)

	def diffuse_reflect(self, surface, cylinder, rng=np.random):
		"""Diffuse reflection: pick a random direction in the hemisphere defined by the surface normal."""
		# Get normal
		if surface == 'side':
			rel = self.position - cylinder.center
			h = np.dot(rel, cylinder.axis)
			proj = cylinder.center + h * cylinder.axis
			n = (self.position - proj)
			n = n / np.linalg.norm(n)
		elif surface == 'top':
			n = cylinder.axis
		elif surface == 'bottom':
			n = -cylinder.axis
		else:
			raise ValueError('Unknown surface for reflection')
		# Sample random direction in hemisphere
		while True:
			v = random_unit_vector(rng)
			if np.dot(v, n) > 0:
				self.direction = v
				break
