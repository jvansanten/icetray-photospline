import glam
import splinetable
import numpy
import Gnuplot
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

print "Reading spline fit..."
coeff = numpy.load(sys.argv[2])
table = splinetable.SplineTable()
table.knots = knots
table.coefficients = coeff
table.periods = periods
table.order = 2

print "Coefficient matrix shape is",coeff.shape

smoothed = glam.grideval(table,axes)
smoothed = numpy.exp(smoothed)-1

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

for y in axes[2]:
	sample = bigdat[bigdat[:,2] == y]
	raw = Gnuplot.Data(sample[:,0],sample[:,1],sample[:,3],title="Raw")
	fit = Gnuplot.Data(sample[:,0],sample[:,1],sample[:,4],title="Fit")
	gp.splot(raw,fit)
	printdiffstats(sample[:,3],sample[:,4])
	print 'NdirC is',y
	raw_input('Press Enter to continue')

