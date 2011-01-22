
import numpy
from .. import photo2numpy, glam, numpy_extensions
import warnings

class Efficiency:
	"""Normalization types from photonics.h"""
	NONE         = 0x00
	RECEIVER     = 0x01 
	SOURCE       = 0x02
	WAVELENGTH   = 0x04
	AREA         = 0x08
	VOLUME       = 0x10
	N_PHOTON     = 0x20
	EMISSION     = 0x40
	USER_DEFINED = 0x80
	DIFFERENTIAL = 0x100

class Geometry:
	SPHERICAL   = 1
	CYLINDRICAL = 2
	CUBIC       = 3

# Class for reading IceCube photonics tables
class photonics_table():
	level = -1

	table       = None
	weights     = None
	bin_centers = None
	bin_widths  = None

	is_integral = False

	filename    = None

	# Constructor. Creates instance and optionally opens pt file.
	def __init__(self, filename=None, normalize=True, mmap=False):
		if filename is not None:
			self.open_file(filename, mmap=mmap)
		if normalize:
			self.normalize()
	    

	# Checks consistency of loaded tables.
	# For now, only checks shapes of various arrays.
	def table_shape_consistent(self):
		shapes = set()

		shapes.add(self.values.shape)
		if self.weights is not None:
			shapes.add(self.weights.shape)
		shapes.add(tuple([len(i) for i in self.bin_centers]))
		shapes.add(tuple([len(i) for i in self.bin_widths]))
		if len(shapes) > 1:
			return 0

		return 1

	# Normalize to absolute bin amplitudes
	def normalize(self):
		eff = self.header['efficiency']
		if not (eff & Efficiency.N_PHOTON):
			# This table still contains raw weights from photomc
			self.values /= self.header['n_photon']
			self.weights /= self.header['n_photon']
			eff = eff | Efficiency.N_PHOTON
		if (eff & Efficiency.DIFFERENTIAL):
			# Someone has made this a dP/dt table. Undo their work.
			if self.values.ndim != 4:
				raise ValueError, "This table is wierd, man."
			shape = [1]*len(self.values.shape)
			shape[-1] = self.values.shape[-1]
			dt = self.bin_widths[-1].reshape(shape)
			self.values *= dt
			eff = eff & ~Efficiency.DIFFERENTIAL

		self.header['efficiency'] = eff

	# Returns number of dimensions in table.
	@property
	def ndim(self):
		return len(self.shape)

	# Returns shape of table.
	@property
	def shape(self):
		if not self.table_shape_consistent():
			raise Exception('Shape consistency check failed')

		return self.values.shape

	def remove_nans_and_infinites(self, dovalues=True, doweights=True):
		if self.weights is not None and doweights:
			self.weights[numpy.logical_not( \
			    numpy.isfinite(self.values))] = 0
		if dovalues:
			self.values [numpy.logical_not( \
			    numpy.isfinite(self.values))] = 0

	@property
	def normed(self):
		"""Has this table been normalized?"""
		if self.values.ndim == 4:
			normval = self.values[:,:,:,-1]
			if (normval[(normval > 0) & \
			    numpy.isfinite(normval)] == 1).all():
				return True
			else:
				return False
		else:
			return True

	def open_file(self, filename, convert=True, mmap=False):
		if mmap:
			warnings.warn("Memory-mapped tables are single-precision. You have been warned.");
			self.filename = filename
		try:
			table = photo2numpy.readl1(filename, mmap)
			self.level = 1
		except ValueError, inst:
			try:
				table = photo2numpy.readl2(filename)
				self.level = 2
			except ValueError, inst:
				print inst
				return 0

		self.values = table[0]

		if table[1] == None:
			self.weights = numpy.ones(self.values.shape)
		else:
			self.weights = table[1]

		self.bin_centers = list(table[2])
		self.bin_widths  = list(table[3])

		self.header = table[4]

		# Level 2 tables get an extra, random element here for
		# some reason
		if len(self.bin_centers) > self.values.ndim:
			self.bin_centers = self.bin_centers[0:self.values.ndim]
		if len(self.bin_widths ) > self.values.ndim:
			self.bin_widths  = self.bin_widths [0:self.values.ndim]
            
		# Check consistency of table shapes and derive type
		ndim = self.ndim
		if ndim == 3:
			self.is_integral = True;

		# Convert to standard format unless user doesn't want this
		if convert:
			self.convert_to_level1()

		return 1

	def convert_to_level1(self):
		if self.level == 0 or self.level == 1:
			return 1

		# For level 2, some axes are reversed.
		if self.level == 2:
			self.values = numpy.rollaxis(self.values, 0, 3)
			if self.weights is not None:
				self.weights = \
				    numpy.rollaxis(self.weights, 0, 3)
			self.bin_centers[2], self.bin_centers[0], \
			    self.bin_centers[1] = self.bin_centers[0], \
			    self.bin_centers[1], self.bin_centers[2]
			self.bin_widths[2],  self.bin_widths[0], \
			    self.bin_widths[1] = self.bin_widths[0], \
			    self.bin_widths[1],  self.bin_widths[2]
            
			from math import pi
			self.bin_centers[1][:] *= 180./pi
			self.bin_widths[1][:]  *= 180./pi
            
			self.level = 1

			return 1
            
		print "Don't know how to convert table with level", self.level
		return 0
            
	def convert_to_level2(self):
		if self.level == 2:
			return 1

		# For level 0/1, some axes are reversed.
		if self.level == 0 or self.level == 1:
			self.values =       numpy.rollaxis(self.values, 2, 0)
			if self.weights is not None:
				self.weights = numpy.rollaxis(self.weights, \
				    2, 0)
			self.bin_centers[0], self.bin_centers[1], \
			    self.bin_centers[2] = self.bin_centers[2], \
			    self.bin_centers[0], self.bin_centers[1]
			self.bin_widths[0], self.bin_widths[1], \
			    self.bin_widths[2] = self.bin_widths[2], \
			    self.bin_widths[0], self.bin_widths[1]

			from math import pi
			self.bin_centers[1][:] *= pi/180.
			self.bin_widths[1][:]  *= pi/180.

			self.level = 2
            
			return 1
            
		print "Don't know how to convert table with level", self.level
		return 0

	def mirror(self,n_rho=0,n_phi=0):
		"""Extend table to rho < 0 and 180 < phi < 360. This may be useful for surpressing edge effects while fitting."""

		if n_rho == 0 and n_phi == 0:
			return None

		if abs(self.bin_widths[1].sum() - 180) > 1e-12:
			raise ValueError, "Only half-cylindrical tables can \
			    be mirrored. Perhaps mirror() has already been \
			    called?"

		## XXX only phi mirroring for now
		new_shape = list(self.values.shape)
		new_shape[0] += n_rho
		new_shape[1] += 2*n_phi

		target_slice = [slice(None)]*self.values.ndim
		source_slice = [slice(None)]*self.values.ndim
		target_slice[0] = slice(n_rho, None)
		target_slice[1] = slice(n_phi, -n_phi)

		# copy values into expanded array
		new_values = numpy.empty(new_shape)
		new_values[target_slice] = self.values

		# replace old values with expanded version
		del self.values
		self.values = new_values

		# copy weights into expanded array
		new_weights = numpy.empty(new_shape)
		new_weights[target_slice] = self.weights

		# replace weights
		del self.weights
		self.weights = new_weights

		# replace bin centers and widths
		for lst in (self.bin_centers, self.bin_widths):
			for i in (0,1):
				new = numpy.empty(new_shape[i])
				new[target_slice[i]] = lst[i]
				lst[i] = new

		# mirror left edge
		source_slice[1] = [2*n_phi - 1 - i for i in xrange(n_phi)]
		target_slice[0] = slice(None)
		target_slice[1] = range(n_phi)
		for array in (self.values, self.weights):
			array[target_slice] = array[source_slice]
		for lst in (self.bin_centers, self.bin_widths):
			lst[1][target_slice[1]] = -(lst[1][source_slice[1]])

		# mirror right edge
		source_slice[1] = [-(2*n_phi - i) for i in xrange(n_phi)]
		target_slice[1] = [-(i+1) for i in xrange(n_phi)]
		for array in (self.values, self.weights):
			array[target_slice] = array[source_slice]
		for lst in (self.bin_centers, self.bin_widths):
			lst[1][target_slice[1]] = 360 - lst[1][source_slice[1]]

		# mirror radial slices
		# negative radii are mirrored, so in reverse order
		source_slice[0] = range(2*n_rho - 1, n_rho - 1, -1)
		target_slice[0] = range(n_rho)
		for lst in (self.bin_centers, self.bin_widths):
			lst[0][target_slice[0]] = -(lst[0][source_slice[0]])

		# mirror the radial slice at each azimuth to negative radii
		for i in xrange(self.bin_centers[1].size):
			# find the opposite slice
			opposite = 180 + self.bin_centers[1][i]
			if opposite > 180: opposite -= 2*(180 - opposite)
			elif opposite < 0: opposite *= -1

			mcenter = abs(opposite - self.bin_centers[1]).argmin()
			source_slice[1] = mcenter
			target_slice[1] = i
			for array in (self.values, self.weights):
				array[target_slice] = array[source_slice]
	
		return None

def melonball(table, weights = None, radius = 1):
	"""Set weights inside a given radius to zero."""
	if table.header['geometry'] != Geometry.CYLINDRICAL:
		raise ValueError, "Can't handle non-cylindrical tables"
	Rho, Z = numpy.meshgrid_nd(table.bin_centers[0], table.bin_centers[2], lex_order=True)
	mask = Rho**2 + Z**2 < radius**2
	if weights is None:
		weights = table.weights
	shape = weights.shape
	for i in xrange(shape[1]):
		if weights.ndim == 3:
			weights[:,i,:][mask] = 0
		else:
			for j in xrange(shape[3]):
				weights[:,i,:,j][mask] = 0
				
def scalp(table, weights = None, low = -820, high = 820):
	"""Set weights outside of the given depth range to zero."""
	
	geo = table.header['geometry']
	depth = table.header['z']
	zenith = table.header['zenith']*numpy.pi/180.0
	if geo == Geometry.CYLINDRICAL:
		Rho, Phi, L = numpy.meshgrid_nd(*table.bin_centers[:3], lex_order=True)
	elif geo == Geometry.SPHERICAL:
		R, Phi, CosPolar = numpy.meshgrid_nd(*table.bin_centers[:3], lex_order=True)
		L = R*CosPolar
		Rho = numpy.sqrt(R**2 - L**2)
		del R, CosPolar
	else:
		raise ValueError, "Unknown geometry type %d" % geo
	Phi *= (numpy.pi/180.0)
		
	z = L*numpy.cos(zenith) + Rho*numpy.sin(zenith)*(numpy.cos(Phi) + numpy.sin(Phi))
	mask = (z > low)&(z < high)
	
	del Rho, Phi, L, z
	
	if weights is None:
		weights = table.weights
	shape = weights.shape
	ddim = len(shape) - len(mask.shape)
	if ddim != 0:
		mask = mask.reshape(mask.shape + (1,)*ddim)
	
	weights[mask] = 0
	