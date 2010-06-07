
import numpy
import photo2numpy

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


# Class for reading IceCube photonics tables
class photonics_table():

    level = -1
    
    table       = None
    weights     = None
    bin_centers = None
    bin_widths  = None

    is_integral = False
    
    filename    = None
    
    # Contructor. Creates instance and optionally opens pt file.
    def __init__(self, filename=None):
        if filename is not None:
            self.open_file(filename)

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
        
    # Returns number of dimensions in table.
    def ndim(self):
        return len(self.shape())
    
    # Returns shape of table.
    def shape(self):
        if not self.table_shape_consistent():
            raise Exception('Shape consistency check failed')
            
        return self.values.shape

    def remove_nans_and_infinites(self, dovalues=True, doweights=True):
        if self.weights is not None \
             and \
           doweights:
            self.weights[numpy.logical_not(numpy.isfinite(self.values))] = 0
        if dovalues:
            self.values [numpy.logical_not(numpy.isfinite(self.values))] = 0

    def open_file(self, filename, convert=True):
        try:
            table = photo2numpy.readl1(filename)
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

        # Level 2 tables get an extra, random element here for some reason
        if len(self.bin_centers) > self.values.ndim:
            self.bin_centers = self.bin_centers[0:self.values.ndim]
        if len(self.bin_widths ) > self.values.ndim:
            self.bin_widths  = self.bin_widths [0:self.values.ndim]
            
        # Check consistency of table shapes and derive type
        ndim = self.ndim()
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
            self.values =       numpy.rollaxis(self.values, 0, 3)
            if self.weights is not None:
                self.weights = numpy.rollaxis(self.weights, 0, 3)
            self.bin_centers[2], self.bin_centers[0], self.bin_centers[1] = \
                self.bin_centers[0], self.bin_centers[1], self.bin_centers[2]
            self.bin_widths[2],  self.bin_widths[0],  self.bin_widths[1] = \
                self.bin_widths[0],  self.bin_widths[1],  self.bin_widths[2]
            
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
                self.weights = numpy.rollaxis(self.weights, 2, 0)
            self.bin_centers[0], self.bin_centers[1], self.bin_centers[2] = \
                self.bin_centers[2], self.bin_centers[0], self.bin_centers[1]
            self.bin_widths[0],  self.bin_widths[1],  self.bin_widths[2] = \
                self.bin_widths[2],  self.bin_widths[0],  self.bin_widths[1]
            
            from math import pi
            self.bin_centers[1][:] *= pi/180.
            self.bin_widths[1][:]  *= pi/180.

            self.level = 2
            
            return 1
            
        print "Don't know how to convert table with level", self.level
        return 0

    def mirror(self,rho=3,phi=3):
	"""Extend table to rho < 0 and 180 < phi < 360. This may be useful for surpressing edge effects while fitting."""
	# XXX only rho mirroring for now
	# the mirror images of the radial slices are in reverse order
	indices = numpy.arange(self.values.shape[1], -1, -1)
	# first rho radial slices in descending radial order	
	mirrorslice = [slice(None,rho,-1),indices] + [slice(None)]*self.values.ndim
	mirror = table.values
	# XXX not done yet...
	return None

    #def __del__(self):
        #
