import numpy
from glam import glam, splinefitstable
import sys

# Hard-coded params

bins = 50
nknots = 15
smooth = 0.1

# Real code

data = numpy.loadtxt(sys.argv[1])
ndim = data.ndim
ranges = numpy.column_stack((data.min(0),data.max(0)))
knots = []
periods = []
print "Axis lengths:"
for r in ranges:
	print "\t",r[0],"-",r[1]
	space = (r[1] - r[0])/nknots
	knots.append(numpy.linspace(r[0]-3.5*space,r[1]+2*space,nknots))
	periods.append(0)

print "Histogramming..."

z,axes = numpy.histogramdd(data,bins=bins,normed=False)
for i in range(0,len(axes)):
	x = axes[i]
	x = x + (x[1] - x[0])/2.
	x.resize(x.size - 1)
	axes[i] = x

print "Loaded histogram with dimensions ",z.shape

print "Beginning spline fit..."
# Set weights equal to the variance equal to the 1 + the observed values
table = glam.fit(numpy.log(z+1.),z + 1.,axes,knots,2,smooth,periods)

print "Saving coefficients to %s..." % (sys.argv[1]+".pspl.fits")
splinefitstable.write(table,sys.argv[1]+".pspl.fits")
