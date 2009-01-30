import numpy

def bspline(knots, x, i, n):
	if n == 0:
		if (x > knots[i] and x < knots[i+1]):
			return 1
		else:
			return 0

	a = (x - knots[i])*bspline(knots, x, i, n-1)/(knots[i+n] - knots[i])
	b = (knots[i+n+1] - x)*bspline(knots, x, i+1, n-1)/(knots[i+n+1] - knots[i+1])
	return a+b

def splinebasis(knots,order,x1):
	splinevals = []
	nsplines = len(knots)-order-1
	i = 0
	while i < nsplines:
		splinevals.append([bspline(knots,x,i,order) for x in x1])
		i = i+1

	return numpy.matrix(numpy.column_stack(splinevals))

