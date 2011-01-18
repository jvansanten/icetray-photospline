from icecube.photospline.glam import glam
from icecube.photospline import splinetable, splinefitstable
import numpy
import Gnuplot
import sys

# Hard-coded params
bins = 50

# Real code

data = numpy.loadtxt(sys.argv[1])
ndim = data.ndim
ranges = numpy.column_stack((data.min(0),data.max(0)))

print "Histogramming..."
z,axes = numpy.histogramdd(data,bins=bins,normed=True)
for i in range(0,len(axes)):
	x = axes[i]
	x = x + (x[1] - x[0])/2.
	x.resize(x.size - 1)
	axes[i] = x

print "Loaded histogram with dimensions ",z.shape

print "Reading spline fit..."
table = splinefitstable.read(sys.argv[2])

print "Coefficient matrix shape is",table.coefficients.shape

smoothed = glam.grideval(table,axes)

# Invert the link function we chose before
smoothed = numpy.exp(smoothed + numpy.log(1e-9))

def printdiffstats(a,b):
	print "Fit Statistics:"
	print "\tMaximum Deviation from Data:",numpy.max(numpy.abs(a - b))
	print "\tRMS Deviation from Data:",numpy.sqrt(numpy.mean((a - b)**2))
	print "\tMax Fractional Deviation from Data:",numpy.max(numpy.abs((a - b)/b))
	print "\tMean Fractional Deviation from Data:",numpy.mean(numpy.abs((a - b)/b))

coords = []
def append_all_perms(a,i):
	global coords
	if i == len(axes)-1:
		for x in axes[i]:
			b = a[:]
			b.append(x)
			coords.append(b)
	else:
		for x in axes[i]:
			b = a[:]
			b.append(x)
			append_all_perms(b,i+1)

append_all_perms([],0)
coord = numpy.asarray(coords)

zvec = z.reshape(z.size)
bigdat = numpy.column_stack((coord,zvec))
bigdat = numpy.column_stack((bigdat,smoothed.reshape(smoothed.size)))

gp = Gnuplot.Gnuplot()
gp.interact()

for y in axes[2]:
	sample = bigdat[bigdat[:,2] == y]
	raw = Gnuplot.Data(sample[:,0],sample[:,1],sample[:,3],title="Raw")
	fit = Gnuplot.Data(sample[:,0],sample[:,1],sample[:,4],title="Fit")
	gp.splot(raw,fit)
	#gp.splot(fit)
	printdiffstats(sample[:,3],sample[:,4])
	print 'NdirC is',y
	raw_input('Press Enter to continue')

