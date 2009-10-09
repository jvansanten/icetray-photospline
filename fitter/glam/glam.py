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

def fit(z,w,coords,knots,order,smooth,periods,bases=None):
	ndim=z.ndim

	table = splinetable.SplineTable()
	table.knots = knots
	table.order = order
	table.periods = periods

	nsplines = []
	for i in range(0,len(knots)):
		if (periods[i] == 0):
			nsplines.append(len(knots[i])-order-1)
		else:
			nsplines.append(len(knots[i]))		

	print "Calculating spline basis..."

	if bases == None:
		Basis = [splinebasis(knots[i],order,coords[i],periods[i]) for i in range(0,ndim)]
	else:
		Basis = [numpy.matrix(i) for i in bases]

	print "Calculating penalty matrix..."

	def calcP(nsplines, dim):
		nspl = nsplines[dim]

		D = numpy.eye(nspl-2,nspl,dtype=float,k=0) + \
		    -2*numpy.eye(nspl-2,nspl,dtype=float,k=1) + \
		    numpy.eye(nspl-2,nspl,dtype=float,k=2)
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

	P = calcP(nsplines,0)
	for i in range(1,ndim):
		P = P + calcP(nsplines,i)
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

def grideval(table, coords):
	results = table.coefficients
	Basis = [splinebasis(table.knots[i], table.order,coords[i],
	    table.periods[i]) for i in range(0,len(table.knots))]

	for i in range(0,results.ndim):
		results = rho(Basis[i].transpose(),results,i)

	return results

