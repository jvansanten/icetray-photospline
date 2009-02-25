from glam import glam, splinetable, splinefitstable
import numpy
import sys
import Gnuplot

# Real code

rawdata = numpy.loadtxt(sys.argv[1])
rawdata = rawdata[rawdata[:,1] == 45]
z = rawdata[:,6]
data = rawdata[:,[0,2,5]]
ranges = numpy.column_stack((data.min(0),data.max(0)))
print "Axis lengths:"
for r in ranges:
	print "\t",r[0],"-",r[1]

munge = [numpy.unique(data[:,i]) for i in range(0,3)]
z = z.reshape(munge[0].size,munge[1].size,munge[2].size)

z = numpy.log(z)
print "Loaded histogram with dimensions ",z.shape

print "Reading spline fit..."
table = splinefitstable.read(sys.argv[2])

print "Coefficient matrix shape is",coeff.shape

smoothed = glam.grideval(table,munge)

def printdiffstats(a,b):
	print "Fit Statistics:"
	print "\tMaximum Deviation from Data:",numpy.max(numpy.abs(a - b))
	print "\tRMS Deviation from Data:",numpy.sqrt(numpy.mean((a - b)**2))
	print "\tMax Fractional Deviation from Data:",numpy.max(numpy.abs((a - b)/b))
	print "\tMean Fractional Deviation from Data:",numpy.mean(numpy.abs((a - b)/b))

coord = rawdata[:,[0,2,5]]
zvec = z.reshape(z.size)
bigdat = numpy.column_stack((coord,zvec))
bigdat = numpy.column_stack((bigdat,smoothed.reshape(smoothed.size)))

gp = Gnuplot.Gnuplot()

for y in munge[2]:
	sample = bigdat[bigdat[:,2] == y]
	raw = Gnuplot.Data(sample[:,0],sample[:,1],sample[:,3],title="Raw")
	fit = Gnuplot.Data(sample[:,0],sample[:,1],sample[:,4],title="Fit")
	gp.splot(raw,fit)
	print 'Displaying slice at z = %f' % y
	printdiffstats(sample[:,3],sample[:,4])
	raw_input('Press Enter to continue')

