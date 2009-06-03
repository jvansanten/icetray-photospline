#from glam import glam, splinefitstable
from glam import splinefitstable
import numpy
import sys
import os

# Hard-coded params

nknots = [13, 7, 17, 20] # [r, phi, z, t]
#nknots = [9, 7, 10, 13] # [r, phi, z, t]
smooth = 0.05

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

	table = photo2numpy.readl1(sys.argv[1])
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

	if sys.argv[1].endswith(".prob"):
		first_slice = [slice(None)]*(len(z.shape)-1) + [slice(0,1)]
		z = numpy.append(z[first_slice],numpy.diff(z,axis=-1),axis=-1) / bin_widths[-1]

except:
	print "Using photo2numpy failed, falling back on text processing...\n"

	data = numpy.loadtxt(sys.argv[1])
	data = data[:,[0,1,2,5,6]]
	z = data[:,4]
	ndim = 4
	try:
		weights = numpy.loadtxt(sys.argv[1] + ".stats")[:,6]
	except:
		weights = numpy.ones(z.shape)
	# Get coordinates
	bin_centers = [numpy.unique(data[:,i]) for i in range(0,4)]

	# Reshape to proper form
	z = z.reshape(bin_centers[0].size,bin_centers[1].size,bin_centers[2].size,bin_centers[3].size)
	weights = weights.reshape(z.shape)

# Compute knot locations 

rknots = numpy.append(numpy.logspace(-3,3,nknots[0]), [1100, 1200, 1300, 1400, 1500])
thetaknots = numpy.append(-1.+numpy.logspace(-3,0,3),numpy.append(numpy.linspace(5,175,nknots[1]),181. - numpy.logspace(0,-3,5)))
#thetaknots = numpy.linspace(-70,300,16)
zknots = numpy.append(numpy.logspace(-2,3,nknots[2]/2), [1100, 1200, 1300, 1400, 1500])
zknots = numpy.append(numpy.append([-1300,-1200,-1100],-1.*numpy.logspace(3,-2,nknots[2]/2)),zknots)
tknots = numpy.append(numpy.append([-1,-0.5,0],numpy.logspace(0,numpy.log10(7000),nknots[3])), [7100, 7150, 7200, 7300, 7400])

periods = [0,0,0,0]
knots = [rknots, thetaknots, zknots, tknots]

# HACK: Move first and last angular bins to 0 and 180

bin_centers[1][0] = 0
bin_centers[1][bin_centers[1].size - 1] = 180

# Convert the input to log-space and drop any NaNs or infinites from the fit
z = numpy.log(z)
w = weights
w[numpy.logical_not(numpy.isfinite(z))] = 0
z[numpy.logical_not(numpy.isfinite(z))] = 0

print "Loaded histogram with dimensions ",z.shape

if len(sys.argv) < 3:
	outputfile = sys.argv[1]+".pspl.fits"
else:
	outputfile = sys.argv[2]

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

