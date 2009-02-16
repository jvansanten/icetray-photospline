import numpy
import glam
import sys

def fithist(data, bins, nknots, smooth, link):
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

	z,axes = numpy.histogramdd(data,bins=bins,normed=True)
	for i in range(0,len(axes)):
		x = axes[i]
		x = x + (x[1] - x[0])/2.
		x.resize(x.size - 1)
		axes[i] = x

	# Get the actual bin counts for weighting the fit
	counts = numpy.histogramdd(data,bins=bins,normed=False)[0]

	print "Loaded histogram with dimensions ",z.shape

	print "Beginning spline fit..."

	# Set weights proportional to the (Poisson) variance: 1 + counts 
	return glam.fit(link(z),counts + 1.,axes,knots,2,smooth,periods)

