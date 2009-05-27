import numpy
from glam.bspline import *
import Gnuplot

numpts = 20
knots=range(-3,30)
order=3

x1 = numpy.linspace(0,25,numpts)

# Pick a random complicated function to interpolate
z = numpy.random.poisson(numpy.cos(x1)*numpy.cos(x1) + (x1 -3.0)*(x1-3.0) + 10)
#z = numpy.cos(x1)*numpy.cos(x1) + (x1 -3.0)*(x1-3.0) + 10
#z = numpy.ones(x1.shape) + 5. + 4*numpy.sin(x1)
#z = numpy.ones(x1.shape) + 5. + x1**2/6.

gp = Gnuplot.Gnuplot()
xfine = numpy.sort(numpy.random.uniform(0,25,size=1000))
rawdat = Gnuplot.Data(x1, z, title = "Data")

# See if we can jump to the answer weighting by z and shifting the knots
baseknots = x1 + (numpy.max(x1)-numpy.min(x1))/(2.0*numpts)*(order-1)
interpknots = []
for i in range (order,0,-1):
	interpknots.append(baseknots[0] - i*(baseknots[1] - baseknots[0]))
interpknots.extend(baseknots)
interpknots.append(interpknots[len(interpknots)-1] + (interpknots[len(interpknots)-1] - interpknots[len(interpknots)-2]))

splinterp = Gnuplot.Data(xfine, [sum([z[n]*bspline(interpknots, x, n, order) for n in range(0,len(interpknots)-order-1)]) for x in xfine], with_="lines", title = "Direct Spline Interpolation Attempt")
knotcoeff = Gnuplot.Data(baseknots, z, title="Knots and Coefficients")

knots = interpknots
# Do an overparameterized least-squares fit for comparison
A = splinebasis(knots,order,x1)
result = numpy.linalg.lstsq(A, z)

# Plot the least-squares result
spline = Gnuplot.Data(xfine, [sum([result[0][n]*bspline(knots, x, n, order) for n in range(0,len(knots)-order-1)]) for x in xfine], with_="lines", title="Least Squares")

#gp.plot(rawdat,splinterp,spline)
#gp.set_range("yrange",(-1,17))
gp.set_range("xrange",(-1,26))
gp.plot(rawdat,splinterp,spline,knotcoeff)

