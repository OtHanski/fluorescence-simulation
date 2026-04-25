import csv
from tap import Tap
from samplecell import SurfaceProperties

class CoatingPresets:
	@staticmethod
	def fully_absorptive():
		"""All photons absorbed, no reflection or conversion."""
		return SurfaceProperties(
			absorption={'blue': 1.0, 'uv': 1.0},
			conversion={'blue': 0.0, 'uv': 0.0},
			specular={'blue': 0.0, 'uv': 0.0},
			diffuse={'blue': 0.0, 'uv': 0.0}
		)

	@staticmethod
	def specular_mirror():
		"""Highly specular reflective mirror, no absorption or conversion."""
		return SurfaceProperties(
			absorption={'blue': 0.0, 'uv': 0.0},
			conversion={'blue': 0.0, 'uv': 0.0},
			specular={'blue': 1.0, 'uv': 1.0},
			diffuse={'blue': 0.0, 'uv': 0.0}
		)

	@staticmethod
	def diffuse_converter():
		"""Diffuse reflective, 50% UV->blue conversion, 50% diffuse reflection for UV, blue is just reflected diffusely."""
		return SurfaceProperties(
			absorption={'blue': 0.0, 'uv': 0.0},
			conversion={'blue': 0.0, 'uv': 0.5},  # 50% of UV photons convert to blue
			specular={'blue': 0.0, 'uv': 0.0},
			diffuse={'blue': 1.0, 'uv': 0.5}      # blue: 100% diffuse, uv: 50% diffuse (rest converts)
		)
from samplecell import SampleCell
import numpy as np

class Preset:
	@staticmethod
	def build_system():
		"""
		Returns a SampleCell with:
		- Bottom: 16mm radius, 30cm long cylinder
		- Top: 150mm radius, 40cm long cylinder, stacked on bottom
		- 6 sensors radially at the edge of the bottom cylinder
		- 9 sensors at the top flange around the axis of the top cylinder
		"""
		cell = SampleCell()
		# Bottom cylinder
		r1 = 0.016  # 16 mm
		h1 = 0.30   # 30 cm
		c1 = np.array([0,0,0])
		cyl1 = cell.add_cylinder(center=c1, radius=r1, height=h1)
		# Top cylinder
		r2 = 0.150  # 150 mm
		h2 = 0.40   # 40 cm
		cyl2 = cell.add_cylinder(center=None, radius=r2, height=h2, connect_to=cyl1)

		# 6 sensors at edge of bottom cylinder (z=0, r=r1, evenly spaced)
		for i in range(6):
			angle = 2 * np.pi * i / 6
			x = r1 * np.cos(angle)
			y = r1 * np.sin(angle)
			center = [x, y, 0]
			normal = [x, y, 0]
			normal = np.array(normal)
			if np.linalg.norm(normal) == 0:
				normal = [0,0,1]
			else:
				normal = normal / np.linalg.norm(normal)
			cell.add_sensor(center=center, normal=normal)

		# 9 sensors at top flange (z=h1+h2, r=r2, evenly spaced)
		z_top = h1 + h2
		for i in range(9):
			angle = 2 * np.pi * i / 9
			x = r2 * np.cos(angle)
			y = r2 * np.sin(angle)
			center = [x, y, z_top]
			normal = [x, y, 0]
			normal = np.array(normal)
			if np.linalg.norm(normal) == 0:
				normal = [0,0,1]
			else:
				normal = normal / np.linalg.norm(normal)
			cell.add_sensor(center=center, normal=normal)

		return cell
import numpy as np
from photon import Photon

class PhotonSource:
	def __init__(self, center, radius, temperature, mass, g=9.81, axis=(0,0,1), wavelength='uv', rng=np.random):
		"""
		center: (x, y, z) tuple for base center of the cloud
		radius: radius of the cloud (adjustable)
		temperature: temperature in Kelvin
		mass: mass of gas particle (kg)
		g: gravitational acceleration (m/s^2)
		axis: direction of cylinder axis (default z)
		wavelength: photon wavelength to emit
		rng: random number generator
		"""
		self.center = np.array(center, dtype=float)
		self.radius = radius
		self.temperature = temperature
		self.mass = mass
		self.g = g
		self.axis = np.array(axis, dtype=float) / np.linalg.norm(axis)
		self.wavelength = wavelength
		self.rng = rng
		# Calculate height from temperature and gravity (Boltzmann distribution)
		k_B = 1.380649e-23  # Boltzmann constant (J/K)
		# Height where exp(-mgh/kT) ~ 0.01 (i.e., 99% of atoms below this height)
		self.height = -k_B * temperature / (mass * g) * np.log(0.01)

	def random_point(self):
		"""Generate a random point inside the cylindrical cloud."""
		# Random height (Boltzmann distribution)
		k_B = 1.380649e-23
		z = self.rng.uniform(0, self.height)
		# For true Boltzmann, sample z with exp(-m g z / kT) weighting
		# Inverse CDF: z = -kT/(mg) * log(1-u), u~Uniform(0,1)
		u = self.rng.uniform(0,1)
		z = -k_B * self.temperature / (self.mass * self.g) * np.log(1-u)
		# Random radius (uniform in area)
		r = self.radius * np.sqrt(self.rng.uniform(0,1))
		theta = self.rng.uniform(0, 2*np.pi)
		x = r * np.cos(theta)
		y = r * np.sin(theta)
		# Place in 3D
		base = self.center
		# Build orthonormal basis
		z_axis = self.axis
		# Find a perpendicular vector
		if abs(z_axis[0]) < 0.9:
			x_axis = np.cross(z_axis, [1,0,0])
		else:
			x_axis = np.cross(z_axis, [0,1,0])
		x_axis /= np.linalg.norm(x_axis)
		y_axis = np.cross(z_axis, x_axis)
		pos = base + x * x_axis + y * y_axis + z * z_axis
		return pos

	def random_direction(self):
		"""Return a random unit vector (isotropic emission)."""
		phi = self.rng.uniform(0, 2*np.pi)
		costheta = self.rng.uniform(-1, 1)
		sintheta = np.sqrt(1 - costheta**2)
		return np.array([sintheta * np.cos(phi), sintheta * np.sin(phi), costheta])

	def emit_photon(self):
		"""Generate a photon at a random position and direction in the source."""
		pos = self.random_point()
		direction = self.random_direction()
		return Photon(pos, direction, self.wavelength)


# --- Simulation CLI ---

# Default particle mass: Rb-87 atom (kg)
_RB87_MASS = 87 * 1.66053906660e-27
# Cloud radius inside the bottom cylinder (m)
_CLOUD_RADIUS = 0.012

class SimArgs(Tap):
	n_photons: int   # Number of UV photons to simulate
	temperature: float  # Temperature of the gas in Kelvin
	output: str = 'results.csv'  # Output CSV file path
	mass: float = _RB87_MASS  # Gas particle mass in kg
	max_bounces: int = 1000  # Safety limit on bounces per photon
	seed: int = None  # Random seed for reproducibility


def run_simulation(args: SimArgs):
	rng = np.random.RandomState(args.seed)

	# Build cell
	cell = Preset.build_system()

	# Photon source: centred at base of bottom cylinder
	source = PhotonSource(
		center=[0, 0, 0],
		radius=_CLOUD_RADIUS,
		temperature=args.temperature,
		mass=args.mass,
		wavelength='uv',
		rng=rng
	)

	records = []
	for i in range(args.n_photons):
		photon = source.emit_photon()
		outcome = 'lost'
		bounces = 0

		while not photon.absorbed and bounces < args.max_bounces:
			event, info = cell.traverse_photon(photon, rng=rng)
			bounces += 1
			if event == 'detected':
				outcome = 'detected'
				break
			elif event == 'absorbed':
				outcome = 'absorbed'
				break
			elif event == 'escaped':
				outcome = 'lost'
				break
			# specular, diffuse, converted: photon continues

		if bounces >= args.max_bounces:
			outcome = 'lost'

		records.append({
			'photon_id': i,
			'outcome': outcome,
			'x': photon.position[0],
			'y': photon.position[1],
			'z': photon.position[2],
			'wavelength': photon.wavelength,
		})

	with open(args.output, 'w', newline='') as f:
		writer = csv.DictWriter(f, fieldnames=['photon_id', 'outcome', 'x', 'y', 'z', 'wavelength'])
		writer.writeheader()
		writer.writerows(records)

	n_det = sum(1 for r in records if r['outcome'] == 'detected')
	n_abs = sum(1 for r in records if r['outcome'] == 'absorbed')
	n_lost = sum(1 for r in records if r['outcome'] == 'lost')
	print(f"Simulated {args.n_photons} photons.")
	print(f"  Detected:  {n_det}  ({100*n_det/args.n_photons:.1f}%)")
	print(f"  Absorbed:  {n_abs}  ({100*n_abs/args.n_photons:.1f}%)")
	print(f"  Lost:      {n_lost}  ({100*n_lost/args.n_photons:.1f}%)")
	print(f"Results written to: {args.output}")


if __name__ == '__main__':
	args = SimArgs().parse_args()
	run_simulation(args)
