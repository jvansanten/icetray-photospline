import numpy
import splinetable
from bspline import *

def box(A,B):
	ea = numpy.ones((1,A.shape[1]),float)
	eb = numpy.ones((1,B.shape[1]),float)
	
	return numpy.matrix(numpy.asarray(numpy.kron(A, eb)) * \
	    numpy.asarray(numpy.kron(ea, B)))

def rho(A,B,p):
	sa = A.shape
	sb = B.shape

	newaxes = range(p,B.ndim) + range(0,p)

	B = B.transpose(newaxes)
	nonp = numpy.prod([B.shape[i] for i in range(1,B.ndim)])
	B = numpy.reshape(B,(B.shape[0],nonp))
	
	C = numpy.asarray(A.transpose() * numpy.matrix(B))
	C = numpy.reshape(C,[A.shape[1]] + [sb[i] for i in newaxes[1:]])

	# Now invert the permutation we did to B on C
	invaxes = newaxes[:]
	for i in range(0,len(newaxes)): invaxes[newaxes[i]] = i
	C = C.transpose(invaxes)

	return C

def fit(z,w,coords,knots,order,smooth,periods=None,penorder=2,bases=None):
	ndim=z.ndim

	table = splinetable.SplineTable()
	table.knots = knots
	table.order = order

	order = numpy.asarray(order,dtype=long)
	if order.size == 1:
		order = order * numpy.ones(len(knots),dtype=long)
	penorder = numpy.asarray(penorder,dtype=long)
	if penorder.size == 1:
		penorder = penorder * numpy.ones(len(knots),dtype=long)
	if periods == None:
		periods = numpy.zeros(len(knots))
	table.periods = periods

	nsplines = []
	for i in range(0,len(knots)):
		if periods[i] == 0:
			nsplines.append(len(knots[i])-order[i]-1)
		else:
			nsplines.append(len(knots[i]))		

	print "Calculating spline basis..."

	if bases == None:
		Basis = [splinebasis(knots[i],order[i],coords[i],periods[i]) for i in range(0,ndim)]
	else:
		Basis = [numpy.matrix(i) for i in bases]

	print "Calculating penalty matrix..."

	def calcP(nsplines, knots, dim, order, porder):
		nspl = nsplines[dim]
		knots = knots[dim]

		D = numpy.zeros((nspl-order,nspl),dtype=float)
		splcents = knots

		def divdiff(x):
			# Calculate divided difference coefficients
			# in order to estimate derivatives.

			if len(x) == 2:
				return numpy.asarray([-1.0,1.0])/(x[1] - x[0])

			return (len(x)-1.0)*(numpy.append(0,divdiff(x[1:])) -
			    numpy.append(divdiff(x[:-1]),0))/(x[-1] - x[0])

		for i in range(0, len(D)):
			D[i][i:i+porder+1] = divdiff(splcents[i:i+porder+1])

		D = numpy.matrix(D)
		DtD = D.transpose() * D

		def prodterm(i):
			if (i == dim):
				return DtD
			else:
				return numpy.eye(nsplines[i],dtype=float)

		a = prodterm(0)
		i = 1
		while i < ndim:
			b = prodterm(i)
			a = numpy.kron(a,b)
			i = i+1

		return a

	P = calcP(nsplines,knots,0,order[0],penorder[0])
	for i in range(1,ndim):
		P = P + calcP(nsplines,knots,i,order[i],penorder[i])
	P = smooth*P

	sidelen = numpy.product(nsplines)
	a = numpy.reshape(numpy.zeros(sidelen,float),nsplines)

	print "Reticulating splines..."

	n = 0
	while n < 1:
		n = n+1

		F = w
		R = w*z
		for i in range(0,ndim):
			print "\tProcessing dimension",i
			F = rho(box(Basis[i],Basis[i]),F,i)
			R = rho(Basis[i],R,i)

		Fshape = []
		for i in range(0,ndim): Fshape.extend([nsplines[i],nsplines[i]])
		F = numpy.reshape(numpy.asarray(F), Fshape)

		# Now transpose F: first the even axes, then the odd
		Fshape = range(0,F.ndim,2) + range(1,F.ndim,2)
		F = F.transpose(Fshape)

		F = numpy.reshape(F,(sidelen,sidelen))
		r = numpy.reshape(R,(sidelen,1))

		F = F + P

		print "Computing iteration %d least squares solution..." % n

		result = numpy.linalg.lstsq(F, r)
		a = numpy.reshape(result[0],nsplines)

	table.coefficients = a
	return table

def grideval(table, coords, bases=None):
	results = table.coefficients
	order = numpy.asarray(table.order,dtype=long)
	if order.size == 1:
		order = order * numpy.ones(len(table.knots),dtype=long)

	if bases == None:
		Basis = [splinebasis(table.knots[i], order[i],coords[i],
		    table.periods[i]) for i in range(0,len(table.knots))]
	else:
		Basis = bases

	for i in range(0,results.ndim):
		results = rho(Basis[i].transpose(),results,i)

	return results

