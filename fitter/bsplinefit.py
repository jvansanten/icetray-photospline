import numpy
from bspline import *
import Gnuplot

numpts = 500
knots=range(-2,30)
order=2

#x1 = numpy.sort(numpy.random.normal(15,4,size=numpts))
x1 = numpy.sort(numpy.random.uniform(0,25,size=numpts))
z = numpy.random.poisson(numpy.cos(x1)*numpy.cos(x1) + (x1 -3.0)*(x1-3.0) + 10)

splinevals = []
i = 0
while i < len(knots)-order-1:
	splinevals.append([bspline(knots,x,i,order) for x in x1])
	i = i+1

A = numpy.column_stack(splinevals)
result = numpy.linalg.lstsq(A, z)

gp = Gnuplot.Gnuplot()
xfine = numpy.sort(numpy.random.uniform(0,25,size=1000))
rawdat = Gnuplot.Data(x1, z)
spline = Gnuplot.Data(xfine, [sum([result[0][n]*bspline(knots, x, n, order) for n in range(0,len(knots)-2-1)]) for x in xfine], with_="lines")
gp.plot(rawdat,spline)

