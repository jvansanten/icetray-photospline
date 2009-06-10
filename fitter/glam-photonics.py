from glam import splinefitstable
from optparse import OptionParser
import numpy
import sys
import os

# Hard-coded params

nknots = [13, 5, 15, 20] # [r, phi, z, t]
#nknots = [9, 7, 10, 13] # [r, phi, z, t]
smooth = 100

# Parse arguments

usage = "usage: %prog [options] table.pt [output.fits]"
optparser = OptionParser(usage=usage)
optparser.add_option("-e", "--epsilon", dest="epsilon", type="float",
		     help="epsilon with which to replace table zeros")
(opts, args) = optparser.parse_args()

# Real code

# Try to get SPGLAM, fall back on Python GLAM
try:
	from glam import spglam as glam
except:
	print "SPGLAM not found, falling back on Python GLAM...\n"
	from glam import glam

# Try to use photo2numpy to read tables directly, fall back on it being text
try:
	import photo2numpy

	table = photo2numpy.readl1(args[0])
	z = table[0]
	if table[1] == None:
		weights = numpy.ones(z.shape)
	else:
		weights = table[1]
	bin_centers = table[2]
	bin_widths  = table[3]
	ndim = z.ndim

	# Now convert to PDF from CDF if we got a .prob table
	# This uses finite differences, and then divides by the bin widths
	#
	# NB: we can only do this with photo2numpy
	# NB: this is a disgusting hack

	if args[0].endswith(".prob"):
		print "Input table ends in .prob -- attempting to recover PDF"
		first_slice = [slice(None)]*(len(z.shape)-1) + [slice(0,1)]
		z = numpy.append(z[first_slice],numpy.diff(z,axis=-1),axis=-1) / bin_widths[-1]

except:
	print "Using photo2numpy failed, falling back on text processing...\n"

	data = numpy.loadtxt(args[0])
	data = data[:,[0,1,2,5,6]]
	z = data[:,4]
	ndim = 4
	try:
		weights = numpy.loadtxt(args[0] + ".stats")[:,6]
	except:
		weights = numpy.ones(z.shape)
	# Get coordinates
	bin_centers = [numpy.unique(data[:,i]) for i in range(0,4)]

	# Reshape to proper form
	z = z.reshape(bin_centers[0].size,bin_centers[1].size,bin_centers[2].size,bin_centers[3].size)
	weights = weights.reshape(z.shape)

# Compute knot locations using a sparsified version of the bin centers as
#    central knots.

radial_extent = 600

coreknots = [bin_centers[i][numpy.unique(numpy.int32(numpy.linspace(0,len(bin_centers[i])-1,nknots[i])))] for i in range(0,ndim)]

# Now append the extra knots off both ends of the axis in order to provide
# full support at the boundaries

rknots = numpy.append(numpy.append([-1, -0.5, -0.1],coreknots[0]),numpy.asarray([100, 200, 300, 400, 500]) + radial_extent)
thetaknots = numpy.append(numpy.append([-1, -0.5, -0.1],coreknots[1]),[180.1,180.2,180.3,180.4,180.5])
zknots = numpy.append(numpy.append([-800, -700, -600],coreknots[2]),[600,700,800,900,1000])

if ndim==4:
	tknots = numpy.append(numpy.append([-1,-0.5,0],coreknots[3]),[7100, 7150, 7200, 7300, 7400])
	tknots = numpy.sort(numpy.append([1, 5, 10],tknots))
	periods = [0,0,0,0]
	knots = [rknots, thetaknots, zknots, tknots]
else:
	periods = [0,0,0]
	knots = [rknots, thetaknots, zknots]

print 'Number of knots used: ',[len(a) for a in knots]

# HACK: Move first and last angular bins to 0 and 180

bin_centers[1][0] = 0
bin_centers[1][bin_centers[1].size - 1] = 180

# Set up weights and data array
w = weights

if opts.epsilon != None:
	w[z == 0] = 0.01
	z = z + opts.epsilon

# Convert the input to log-space and drop any NaNs or infinites from the fit
z = numpy.log(z)
w[numpy.logical_not(numpy.isfinite(z))] = 0
z[numpy.logical_not(numpy.isfinite(z))] = 0

print "Loaded histogram with dimensions ",z.shape

if len(args) < 3:
	outputfile = args[0]+".pspl.fits"
else:
	outputfile = args[1]

if os.path.exists(outputfile):
	if raw_input("File %s exists. Overwrite? (y/n)" % outputfile) == 'y':
		os.unlink(outputfile)
	else:
		sys.exit()

print "Beginning spline fit..."
table = glam.fit(z,w,bin_centers,knots,2,smooth,periods)

print "Saving table to %s..." % outputfile
splinefitstable.write(table,outputfile)

smoothed = glam.grideval(table,bin_centers)
resid = (smoothed - z)[w != 0]
fracresid = ((smoothed - z)/z)[w != 0]

print "Fit Statistics:"
print "\tMaximum Deviation from Data:",numpy.max(numpy.abs(resid))
print "\tRMS Deviation from Data:",numpy.sqrt(numpy.mean(resid**2))
print "\tMax Fractional Deviation from Data:",numpy.max(numpy.abs(fracresid))
print "\tMean Fractional Deviation from Data:",numpy.mean(numpy.abs(fracresid))

