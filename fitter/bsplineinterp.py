import numpy
from glam.bspline import *
import Gnuplot

numpts = 25
knots=range(-3,30)
order=3

x1 = numpy.linspace(0,25,numpts)

# Pick a random complicated function to interpolate
#z = numpy.random.poisson(numpy.cos(x1)*numpy.cos(x1) + (x1 -3.0)*(x1-3.0) + 10)
z = numpy.cos(x1)*numpy.cos(x1) + (x1 -3.0)*(x1-3.0) + 10

# Do an overparameterized least-squares fit for comparison
A = splinebasis(knots,order,x1)
result = numpy.linalg.lstsq(A, z)

gp = Gnuplot.Gnuplot()
xfine = numpy.sort(numpy.random.uniform(0,25,size=1000))
rawdat = Gnuplot.Data(x1, z)

# Plot the least-squares result
spline = Gnuplot.Data(xfine, [sum([result[0][n]*bspline(knots, x, n, order) for n in range(0,len(knots)-order-1)]) for x in xfine], with_="lines")

# Now see if we can jump to the answer weighting by z and shifting the knots
splinterp = Gnuplot.Data(xfine, [sum([z[n]*bspline(x1 + 0.5*(order-1), x, n-order, order) for n in range(order,len(x1)-1)]) for x in xfine], with_="lines")
gp.plot(rawdat,splinterp,spline)

