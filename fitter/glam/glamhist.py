import numpy
import sys

try:
	import spglam as glam
except:
	print "SPGLAM not found, falling back on Python GLAM..."
	import glam

def fithist(data, weights, bins, nknots, smooth, link, ranges = None, penorder=2, verbose=True):
	ndim = data.ndim
	knots = []
	if ranges == None:
		ranges = numpy.column_stack((data.min(0),data.max(0)))

	nknots = numpy.asarray(nknots,dtype=long)
	if nknots.size == 1:
		nknots = nknots * numpy.ones(len(ranges),dtype=long)

	print "Axis lengths:"
	for i in range(0,len(ranges)):
		r = ranges[i]
		print "\t",r[0],"-",r[1]
		space = (r[1] - r[0])/nknots[i]
		knots.append(numpy.linspace(r[0]-3*space,r[1]+2*space,
		    nknots[i]+5))

	print "Histogramming..."

	z,axes = numpy.histogramdd(data,bins=bins,normed=True,weights=weights)
	for i in range(0,len(axes)):
		x = axes[i]
		x = x + (x[1] - x[0])/2.
		x.resize(x.size - 1)
		axes[i] = x

	# Get the actual bin counts for weighting the fit
	counts = numpy.histogramdd(data,bins=bins,normed=False)[0]

	print "Loaded histogram with dimensions ",z.shape

	# Compute weights and transform data according to the link function
	# Set weights proportional to the (Poisson) variance: 1 + counts 

	z = link(z)
	w = counts + 1e-4

	# Hose data points where the link function blew up, setting their weights to 0 
	w[numpy.isinf(z)] = 0
	w[numpy.isnan(z)] = 0
	z[numpy.isinf(z)] = 0
	z[numpy.isnan(z)] = 0

	print "Beginning spline fit..."

	table = glam.fit(z,w,axes,knots,2,smooth,penorder=penorder)

	if verbose:
		smoothed = glam.grideval(table,axes)
		resid = (smoothed - z)[w > 2]
		fracresid = ((smoothed - z)/z)[w > 2]

		print "Fit Statistics:"
		print "\tMaximum Deviation from Data:",numpy.max(numpy.abs(resid))
		print "\tRMS Deviation from Data:",numpy.sqrt(numpy.mean(resid**2))
		print "\tMax Fractional Deviation from Data:",numpy.max(numpy.abs(fracresid))
		print "\tMean Fractional Deviation from Data:",numpy.mean(numpy.abs(fracresid))

	return table

