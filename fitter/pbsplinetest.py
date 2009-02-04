import numpy
import Gnuplot
import bspline

# Test of periodic bspline basis generator

knots = numpy.linspace(2,10,5)
x = numpy.random.uniform(0,10,size=200)
x = numpy.sort(x)

gp = Gnuplot.Gnuplot()
spline = [Gnuplot.Data(x, [bspline.pbspline(knots, x_ele, n, 2,10) for x_ele in x], with_="lines") for n in range(0,len(knots))]
gp.plot(spline[0],spline[1],spline[2],spline[3])
