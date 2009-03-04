from glam import glam, splinefitstable
import numpy
import sys
import os

# Hard-coded params

nknots = 15
smooth = 0.01

# Real code

rawdata = numpy.loadtxt(sys.argv[1])
#rawdata = rawdata[rawdata[:,1] == 45]
data = rawdata[:,[0,1,2,6]]
ndim = 3
ranges = numpy.column_stack((data.min(0),data.max(0)))

# Compute knot locations 

rknots = numpy.append(numpy.logspace(-3,3,nknots), [1100, 1200, 1300, 1400, 1500])
thetaknots = numpy.append(-1.+numpy.logspace(-3,0,3),numpy.append(numpy.linspace(5,175,nknots),181. - numpy.logspace(0,-3,5)))
#thetaknots = numpy.linspace(-70,300,16)
zknots = numpy.append(numpy.logspace(-2,3,nknots/2), [1100, 1200, 1300, 1400, 1500])
zknots = numpy.append(numpy.append([-1300,-1200,-1100],-1.*numpy.logspace(3,-2,nknots/2)),zknots)

periods = [0,0,0]
knots = [rknots, thetaknots, zknots]

# Get coordinates

z = data[:,3]
munge = [numpy.unique(data[:,i]) for i in range(0,3)]

# HACK: Move first and last angular bins to 0 and 180

munge[1][0] = 0
munge[1][munge[1].size - 1] = 180

z = z.reshape(munge[0].size,munge[1].size,munge[2].size)

z = numpy.log(z)
w = numpy.ones(z.shape)
w[numpy.isinf(z)] = 0
z[numpy.isinf(z)] = 0

print "Loaded histogram with dimensions ",z.shape

outputfile = sys.argv[1]+".pspl.fits";
if os.path.exists(outputfile):
	if raw_input("File %s exists. Overwrite? (y/n)" % outputfile) == 'y':
		os.unlink(outputfile)
	else:
		sys.exit()

print "Beginning spline fit..."
table = glam.fit(z,w,munge,knots,2,smooth,periods)

print "Saving table to %s..." % outputfile
splinefitstable.write(table,outputfile)

smoothed = glam.grideval(table,munge)

print "Fit Statistics:"
print "\tMaximum Deviation from Data:",numpy.max(numpy.abs(smoothed - z))
print "\tRMS Deviation from Data:",numpy.sqrt(numpy.mean((smoothed - z)**2))
print "\tMax Fractional Deviation from Data:",numpy.max(numpy.abs((smoothed - z)/z))
print "\tMean Fractional Deviation from Data:",numpy.mean(numpy.abs((smoothed - z)/z))
