#!/usr/bin/env python

import numpy
from icecube.photospline.splinetable import SplineTable
from icecube.photospline.splinefitstable import write

spline = SplineTable()
spline.ndim = 2
spline.order = [2 for i in range(spline.ndim)]
spline.knots = [numpy.linspace(0, 1, 10) for i in range(spline.ndim)]
nsplines = tuple(knots.size-order-1 for knots,order in zip(spline.knots, spline.order))
spline.coefficients = numpy.ones(nsplines)

write(spline, "constant.fits")
