import glam
import numpy
import sys

# Hard-coded params

bins = 15
smooth = 0.1

# Real code

data = numpy.loadtxt(sys.argv[1])
ndim = data.ndim
ranges = numpy.column_stack((data.min(0),data.max(0)))
print "Axis lengths:"
for r in ranges:
	print "\t",r[0],"-",r[1]

print "Histogramming..."

z,axes = numpy.histogramdd(data,bins=bins,normed=False)
for i in range(0,len(axes)):
	x = axes[i]
	x = x + (x[1] - x[0])/2.
	x.resize(x.size - 1)
	axes[i] = x

print "Loaded histogram with dimensions ",z.shape

print "Beginning spline fit..."
coeff = glam.fit(z,axes,axes,2,smooth)

print "Saving coefficients to %s..." % (sys.argv[1]+".pspl")
numpy.save(sys.argv[1]+".pspl",coeff)
