import glam
import numpy
import sys

# Hard-coded params

nknots = 25
smooth = 0.1

# Real code

rawdata = numpy.loadtxt(sys.argv[1])
z = rawdata[:,6]
data = rawdata[:,[0,1,2]]
ndim = 3
ranges = numpy.column_stack((data.min(0),data.max(0)))
knots = []
periods = []
print "Axis lengths:"
for r in ranges:
	print "\t",r[0],"-",r[1]
	space = (r[1] - r[0])/nknots
	knots.append(numpy.linspace(r[0]-3.5*space,r[1]+2*space,nknots))
	periods.append(0)

munge = [numpy.unique(data[:,i]) for i in range(0,3)]
z = z.reshape(munge[0].size,munge[1].size,munge[2].size)

z = numpy.log(z)
w = numpy.ones(z.shape)

print "Loaded histogram with dimensions ",z.shape

print "Beginning spline fit..."
table = glam.fit(z,w,munge,knots,2,smooth,periods)
coeff = table.coefficients

print "Saving coefficients to %s..." % (sys.argv[1]+".pspl")
numpy.savetxt(sys.argv[1]+".pspl",numpy.reshape(coeff,coeff.size))
print "Coefficient matrix shape is",coeff.shape

smoothed = glam.grideval(table,munge)

print "Fit Statistics:"
print "\tMaximum Deviation from Data:",numpy.max(numpy.abs(smoothed - z))
print "\tRMS Deviation from Data:",numpy.sqrt(numpy.mean((smoothed - z)**2))
print "\tMax Fractional Deviation from Data:",numpy.max(numpy.abs((smoothed - z)/z))
print "\tMean Fractional Deviation from Data:",numpy.mean(numpy.abs((smoothed - z)/z))
