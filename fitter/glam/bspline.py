import numpy

def pbspline(knots, x, i, n, period):
	if n == 0:
		diff = (x - knots[i % len(knots)]) % period
		width = (knots[(i+1) % len(knots)] - knots[i % len(knots)]) % period
		if diff < width:
			return 1
		else:
			return 0

	span = (knots[(i+n) % len(knots)] - knots[i % len(knots)]) % period
	diff = (x - knots[i % len(knots)]) % period

	a = diff*pbspline(knots, x, i, n-1, period)/span

	span = (knots[(i+n+1) % len(knots)] - knots[(i+1) % len(knots)]) % period
	diff = (knots[(i+n+1) % len(knots)] - x) % period

	b = diff*pbspline(knots, x, i+1, n-1, period)/span
	return a+b

def bspline(knots, x, i, n):
	if n == 0:
		if (x >= knots[i] and x < knots[i+1]):
			return 1
		else:
			return 0

	a = (x - knots[i])*bspline(knots, x, i, n-1)/(knots[i+n] - knots[i])
	b = (knots[i+n+1] - x)*bspline(knots, x, i+1, n-1)/(knots[i+n+1] - knots[i+1])
	return a+b

def splinebasis(knots,order,x1,period = 0):
	splinevals = []
	if period == 0:
		nsplines = len(knots)-order-1
	else:
		nsplines = len(knots)
	i = 0
	while i < nsplines:
		if period == 0:
			splinevals.append([bspline(knots,x,i,order) for x in x1])
		else:
			splinevals.append([pbspline(knots,x,i,order,period) for x in x1])
		i = i+1

	return numpy.matrix(numpy.column_stack(splinevals))

