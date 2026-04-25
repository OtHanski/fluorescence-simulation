import numpy as np

import numpy as np

# --- Sensor class ---
class Sensor:
    def __init__(self, center, normal, width=0.005, height=0.005, efficiency=None):
        self.center = np.array(center, dtype=float)
        self.normal = np.array(normal, dtype=float) / np.linalg.norm(normal)
        self.width = width
        self.height = height
        self.efficiency = efficiency or {'blue': 1.0, 'uv': 1.0}
        if abs(self.normal[0]) < 0.9:
            x_axis = np.cross(self.normal, [1,0,0])
        else:
            x_axis = np.cross(self.normal, [0,1,0])
        x_axis /= np.linalg.norm(x_axis)
        y_axis = np.cross(self.normal, x_axis)
        self.x_axis = x_axis
        self.y_axis = y_axis

    def contains(self, point):
        rel = np.array(point, dtype=float) - self.center
        x = np.dot(rel, self.x_axis)
        y = np.dot(rel, self.y_axis)
        return abs(x) <= self.width/2 and abs(y) <= self.height/2

    def detect(self, photon, rng=np.random):
        eff = self.efficiency.get(photon.wavelength, 0.0)
        return rng.rand() < eff

# --- SurfaceProperties class ---
class SurfaceProperties:
    def __init__(self, absorption, conversion, specular, diffuse):
        self.absorption = absorption
        self.conversion = conversion
        self.specular = specular
        self.diffuse = diffuse

    def get(self, prop, wavelength, coord):
        val = getattr(self, prop)[wavelength]
        if callable(val):
            return val(coord)
        return val

# --- Cylinder class ---
class Cylinder:
    def __init__(self, center, radius, height, axis=(0, 0, 1), connect_to=None,
                 side_properties=None, top_properties=None, bottom_properties=None):
        self.center = np.array(center, dtype=float)
        self.radius = radius
        self.height = height
        self.axis = np.array(axis, dtype=float) / np.linalg.norm(axis)
        self.connect_to = connect_to
        self.side_properties = side_properties or SurfaceProperties(
            {'blue': 0.0, 'uv': 0.0}, {'blue': 0.0, 'uv': 0.0}, {'blue': 1.0, 'uv': 1.0}, {'blue': 0.0, 'uv': 0.0})
        self.top_properties = top_properties or SurfaceProperties(
            {'blue': 0.0, 'uv': 0.0}, {'blue': 0.0, 'uv': 0.0}, {'blue': 1.0, 'uv': 1.0}, {'blue': 0.0, 'uv': 0.0})
        self.bottom_properties = bottom_properties or SurfaceProperties(
            {'blue': 0.0, 'uv': 0.0}, {'blue': 0.0, 'uv': 0.0}, {'blue': 1.0, 'uv': 1.0}, {'blue': 0.0, 'uv': 0.0})

    def get_side_property(self, prop, wavelength, z):
        return self.side_properties.get(prop, wavelength, z)

    def get_endcap_property(self, which, prop, wavelength, r):
        if which == 'top':
            return self.top_properties.get(prop, wavelength, r)
        elif which == 'bottom':
            return self.bottom_properties.get(prop, wavelength, r)
        else:
            raise ValueError("which must be 'top' or 'bottom'")

    def top_center(self):
        return self.center + self.axis * self.height

    def contains(self, point):
        point = np.array(point, dtype=float)
        rel = point - self.center
        h = np.dot(rel, self.axis)
        if 0 <= h <= self.height:
            proj = self.center + h * self.axis
            r = np.linalg.norm(point - proj)
            return r <= self.radius
        return False

# --- SampleCell class ---
class SampleCell:
    def __init__(self):
        self.cylinders = []
        self.connections = []

    def add_cylinder(self, center, radius, height, axis=(0, 0, 1), connect_to=None,
                     side_properties=None, top_properties=None, bottom_properties=None):
        if connect_to is not None:
            center = np.array(connect_to.top_center())
        cyl = Cylinder(center, radius, height, axis, connect_to,
                       side_properties, top_properties, bottom_properties)
        self.cylinders.append(cyl)
        if connect_to is not None:
            self.connections.append((connect_to, cyl))
        return cyl

    def add_sensor(self, center, normal, width=0.005, height=0.005, efficiency=None):
        if not hasattr(self, 'sensors'):
            self.sensors = []
        sensor = Sensor(center, normal, width, height, efficiency)
        self.sensors.append(sensor)
        return sensor

    def check_sensors(self, photon, rng=np.random):
        if not hasattr(self, 'sensors'):
            return None, False
        for sensor in self.sensors:
            if sensor.contains(photon.position):
                detected = sensor.detect(photon, rng)
                return sensor, detected
        return None, False

    def traverse_photon(self, photon, rng=np.random):
        # ...existing code for traverse_photon...
        # (leave as is)

    def contains(self, point):
        return any(cyl.contains(point) for cyl in self.cylinders)

    def get_accessible_cylinders(self, start_cyl=None):
        if not self.cylinders:
            return []
        if start_cyl is None:
            start_cyl = self.cylinders[0]
        visited = set()
        stack = [start_cyl]
        while stack:
            cyl = stack.pop()
            if cyl in visited:
                continue
            visited.add(cyl)
            for (from_cyl, to_cyl) in self.connections:
                if from_cyl == cyl and to_cyl not in visited:
                    stack.append(to_cyl)
        return visited

    def is_accessible(self, point):
        accessible_cyls = self.get_accessible_cylinders()
        return any(cyl.contains(point) for cyl in accessible_cyls)
class SurfaceProperties:
        def __init__(self, absorption, conversion, specular, diffuse):
                """
                Each property can be:
                    - a dict with keys 'blue', 'uv' and values as floats (flat value)
                    - a dict with keys 'blue', 'uv' and values as callables (functions)
                The callables should have the following signatures:
                    - For side: f(z) where z is position along the cylinder axis (0=base, height=top)
                    - For end caps: f(r) where r is the distance from center (0=axis, radius=outer edge)
                """
                self.absorption = absorption
                self.conversion = conversion
                self.specular = specular
                self.diffuse = diffuse

        def get(self, prop, wavelength, coord):
                """
                Evaluate the property (absorption, etc.) for a given wavelength and coordinate.
                prop: 'absorption', 'conversion', 'specular', 'diffuse'
                wavelength: 'blue' or 'uv'
                coord: z (for side) or r (for end cap)
                """
                val = getattr(self, prop)[wavelength]
                if callable(val):
                        return val(coord)
                return val

# import numpy for vector math
import numpy as np

class Cylinder:
    def __init__(self, center, radius, height, axis=(0, 0, 1), connect_to=None,
                 side_properties=None, top_properties=None, bottom_properties=None):
        """
        center: (x, y, z) tuple for the center of the cylinder base
        radius: radius of the cylinder
        height: height of the cylinder
        axis: direction vector of the cylinder axis (default is z-axis)
        connect_to: Cylinder this one is connected to (None if base)
        side_properties, top_properties, bottom_properties: SurfaceProperties for each surface
        """
        self.center = np.array(center, dtype=float)
        self.radius = radius
        self.height = height
        self.axis = np.array(axis, dtype=float) / np.linalg.norm(axis)
        self.connect_to = connect_to  # Reference to another Cylinder
        # Surface properties for each surface
        self.side_properties = side_properties or SurfaceProperties(
            {'blue': 0.0, 'uv': 0.0}, {'blue': 0.0, 'uv': 0.0}, {'blue': 1.0, 'uv': 1.0}, {'blue': 0.0, 'uv': 0.0})
        self.top_properties = top_properties or SurfaceProperties(
            {'blue': 0.0, 'uv': 0.0}, {'blue': 0.0, 'uv': 0.0}, {'blue': 1.0, 'uv': 1.0}, {'blue': 0.0, 'uv': 0.0})
        self.bottom_properties = bottom_properties or SurfaceProperties(
            {'blue': 0.0, 'uv': 0.0}, {'blue': 0.0, 'uv': 0.0}, {'blue': 1.0, 'uv': 1.0}, {'blue': 0.0, 'uv': 0.0})

    def get_side_property(self, prop, wavelength, z):
        """Evaluate a side property at position z along the axis (0=base, height=top)."""
        return self.side_properties.get(prop, wavelength, z)

    def get_endcap_property(self, which, prop, wavelength, r):
        """
        Evaluate an end cap property at radius r (0=center, radius=outer edge).
        which: 'top' or 'bottom'
        """
        if which == 'top':
            return self.top_properties.get(prop, wavelength, r)
        elif which == 'bottom':
            return self.bottom_properties.get(prop, wavelength, r)
        else:
            raise ValueError("which must be 'top' or 'bottom'")

    def top_center(self):
        """Return the center of the top base of the cylinder."""
        return self.center + self.axis * self.height

    def contains(self, point):
        """Check if a point is inside the cylinder."""
        point = np.array(point, dtype=float)
        rel = point - self.center
        h = np.dot(rel, self.axis)
        if 0 <= h <= self.height:
            proj = self.center + h * self.axis
            r = np.linalg.norm(point - proj)
            return r <= self.radius
        return False

class SampleCell:
    def add_sensor(self, center, normal, width=0.005, height=0.005, efficiency=None):
        """Add a sensor to the sample cell."""
        if not hasattr(self, 'sensors'):
            self.sensors = []
        sensor = Sensor(center, normal, width, height, efficiency)
        self.sensors.append(sensor)
        return sensor

    def check_sensors(self, photon, rng=np.random):
        """
        Check if the photon hits any sensor and is detected.
        Returns (sensor, detected) or (None, False).
        """
        if not hasattr(self, 'sensors'):
            return None, False
        for sensor in self.sensors:
            if sensor.contains(photon.position):
                detected = sensor.detect(photon, rng)
                return sensor, detected
        return None, False
    def traverse_photon(self, photon, rng=np.random):
        """
        Move the photon until it interacts with a surface, then simulate the interaction.
        Returns a tuple (event, info) where event is one of 'absorbed', 'specular', 'diffuse', 'converted', or 'escaped'.
        info contains details about the event.
        """
        # Find first collision
        start = photon.position
        direction = photon.direction
        min_dist = np.inf
        hit_cyl = None
        hit_surface = None
        hit_point = None
        for cyl in self.cylinders:
            # Side intersection
            d = direction
            p = start
            ca = cyl.axis
            oc = p - cyl.center
            ca_dot_d = np.dot(ca, d)
            ca_dot_oc = np.dot(ca, oc)
            d_proj = d - ca * ca_dot_d
            oc_proj = oc - ca * ca_dot_oc
            A = np.dot(d_proj, d_proj)
            B = 2 * np.dot(d_proj, oc_proj)
            C = np.dot(oc_proj, oc_proj) - cyl.radius ** 2
            disc = B * B - 4 * A * C
            if A != 0 and disc >= 0:
                sqrt_disc = np.sqrt(disc)
                t1 = (-B - sqrt_disc) / (2 * A)
                t2 = (-B + sqrt_disc) / (2 * A)
                for t in [t1, t2]:
                    if t < 1e-8:
                        continue
                    pt = p + t * d
                    rel = pt - cyl.center
                    h = np.dot(rel, ca)
                    if 0 <= h <= cyl.height:
                        if t < min_dist:
                            min_dist = t
                            hit_cyl = cyl
                            hit_surface = 'side'
                            hit_point = pt
            # End cap intersections (bottom and top)
            for which, cap_center, cap_sign in [
                ('bottom', cyl.center, 0),
                ('top', cyl.top_center(), cyl.height)
            ]:
                denom = np.dot(d, ca)
                if abs(denom) < 1e-8:
                    continue
                t = np.dot(cap_center - p, ca) / denom
                if t < 1e-8:
                    continue
                pt = p + t * d
                rel = pt - cyl.center
                h = np.dot(rel, ca)
                if which == 'bottom' and abs(h) > 1e-8:
                    continue
                if which == 'top' and abs(h - cyl.height) > 1e-8:
                    continue
                # Check if within radius
                proj = pt - cyl.center - ca * h
                r = np.linalg.norm(proj)
                if r <= cyl.radius:
                    if t < min_dist:
                        min_dist = t
                        hit_cyl = cyl
                        hit_surface = which
                        hit_point = pt

        if hit_cyl is None:
            return 'escaped', {'position': photon.position.copy()}

        # Move photon to collision point
        photon.position = hit_point

        # Check for sensor hit
        sensor, detected = self.check_sensors(photon, rng)
        if sensor is not None and detected:
            photon.absorb()
            return 'detected', {'position': photon.position.copy(), 'sensor': sensor}

        # Get surface properties
        if hit_surface == 'side':
            z = np.dot(hit_point - hit_cyl.center, hit_cyl.axis)
            props = hit_cyl.side_properties
            coord = z
        else:
            rel = hit_point - hit_cyl.center
            h = np.dot(rel, hit_cyl.axis)
            r = np.linalg.norm(rel - hit_cyl.axis * h)
            props = hit_cyl.top_properties if hit_surface == 'top' else hit_cyl.bottom_properties
            coord = r
        wl = photon.wavelength
        # Probabilities
        p_abs = props.get('absorption', wl, coord)
        p_conv = props.get('conversion', wl, coord)
        p_spec = props.get('specular', wl, coord)
        p_diff = props.get('diffuse', wl, coord)
        # Normalize reflection probabilities
        p_refl = p_spec + p_diff
        if p_refl > 1.0:
            p_spec /= p_refl
            p_diff /= p_refl
            p_refl = 1.0
        # Decide event
        x = rng.rand()
        if x < p_abs:
            photon.absorb()
            return 'absorbed', {'position': photon.position.copy(), 'surface': hit_surface}
        x -= p_abs
        if x < p_conv:
            photon.convert_wavelength()
            return 'converted', {'position': photon.position.copy(), 'surface': hit_surface, 'new_wavelength': photon.wavelength}
        x -= p_conv
        if x < p_spec:
            photon.specular_reflect(hit_surface, hit_cyl)
            return 'specular', {'position': photon.position.copy(), 'surface': hit_surface, 'direction': photon.direction.copy()}
        x -= p_spec
        if x < p_diff:
            photon.diffuse_reflect(hit_surface, hit_cyl, rng)
            return 'diffuse', {'position': photon.position.copy(), 'surface': hit_surface, 'direction': photon.direction.copy()}
        # If nothing else, photon escapes
        return 'escaped', {'position': photon.position.copy()}
    def __init__(self):
        self.cylinders = []
        self.connections = []  # List of (from_cyl, to_cyl) tuples

    def add_cylinder(self, center, radius, height, axis=(0, 0, 1), connect_to=None,
                     side_properties=None, top_properties=None, bottom_properties=None):
        """
        Add a cylinder. If connect_to is given, connect this cylinder's base to the top of connect_to.
        Surface properties can be specified for side, top, and bottom.
        """
        if connect_to is not None:
            # Place base at the top of connect_to
            center = np.array(connect_to.top_center())
        cyl = Cylinder(center, radius, height, axis, connect_to,
                       side_properties, top_properties, bottom_properties)
        self.cylinders.append(cyl)
        if connect_to is not None:
            self.connections.append((connect_to, cyl))
        return cyl

    def contains(self, point):
        """
        Check if a point is inside any of the cylinders (closed system).
        """
        return any(cyl.contains(point) for cyl in self.cylinders)

    def get_accessible_cylinders(self, start_cyl=None):
        """
        Return all cylinders accessible from the base (connected component).
        """
        if not self.cylinders:
            return []
        # If no start_cyl, use the first (base) cylinder
        if start_cyl is None:
            start_cyl = self.cylinders[0]
        visited = set()
        stack = [start_cyl]
        while stack:
            cyl = stack.pop()
            if cyl in visited:
                continue
            visited.add(cyl)
            # Add all directly connected cylinders
            for (from_cyl, to_cyl) in self.connections:
                if from_cyl == cyl and to_cyl not in visited:
                    stack.append(to_cyl)
        return visited

    def is_accessible(self, point):
        """
        Check if a point is inside any cylinder and accessible from the base.
        """
        accessible_cyls = self.get_accessible_cylinders()
        return any(cyl.contains(point) for cyl in accessible_cyls)
import numpy as np

class Cylinder:
    def __init__(self, center, radius, height, axis=(0, 0, 1)):
        """
        center: (x, y, z) tuple for the center of the cylinder base
        radius: radius of the cylinder
        height: height of the cylinder
        axis: direction vector of the cylinder axis (default is z-axis)
        """
        self.center = np.array(center)
        self.radius = radius
        self.height = height
        self.axis = np.array(axis) / np.linalg.norm(axis)

    def contains(self, point):
        """
        Check if a point is inside the cylinder.
        """
        point = np.array(point)
        # Project point onto axis
        rel = point - self.center
        h = np.dot(rel, self.axis)
        if 0 <= h <= self.height:
            # Distance from axis
            proj = self.center + h * self.axis
            r = np.linalg.norm(point - proj)
            return r <= self.radius
        return False

class SampleCell:
    def __init__(self):
        self.cylinders = []

    def add_cylinder(self, center, radius, height, axis=(0, 0, 1)):
        cyl = Cylinder(center, radius, height, axis)
        self.cylinders.append(cyl)

    def contains(self, point):
        """
        Check if a point is inside any of the cylinders in the sample cell.
        """
        return any(cyl.contains(point) for cyl in self.cylinders)