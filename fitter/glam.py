import numpy
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


def fit(z,coords,knots,order,smooth):
	ndim=z.ndim

	nsplines = [len(knots[i])-order-1 for i in range(0,ndim)]

	def box(A,B):
		ea = numpy.ones((1,A.shape[1]),float)
		eb = numpy.ones((1,B.shape[1]),float)
		
		return numpy.matrix(numpy.asarray(numpy.kron(A, eb)) * \
		    numpy.asarray(numpy.kron(ea, B)))

	print "Calculating spline basis..."

	Basis = [splinebasis(knots[i],order,coords[i]) for i in range(0,ndim)]

	w = z
	#w = numpy.ones(z.shape)
	#w = 1.0/z
	#w[z == 0.0] = 0.0005

	print "Calculating penalty matrix..."

	def calcP(nsplines, dim):
		D = numpy.matrix(numpy.diag(numpy.ones(nsplines[dim],float),0) + \
		    numpy.diag(-2*numpy.ones(nsplines[dim]-1,float),-1) + \
		    numpy.diag(numpy.ones(nsplines[dim]-2,float),-2))
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
		P = P + calcP(nsplines,1)
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

		print "Computing iteration %d least squares solution..." % n

		result = numpy.linalg.lstsq(F + P, r)
		a = numpy.reshape(result[0],nsplines)

	return a

def smootheddata(coeffs, knots, order, coords):
	results = coeffs
	Basis = [splinebasis(knots[i],order,coords[i]) for i in range(0,len(knots))]
	for i in range(0,results.ndim):
		results = rho(Basis[i].transpose(),results,i)

	return results

