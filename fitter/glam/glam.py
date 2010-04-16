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

def fit(z,w,coords,knots,order,smooth,periods=None,penalties={2:None},bases=None,iterations=1):
	ndim=z.ndim

	table = splinetable.SplineTable()
	table.knots = knots
	table.order = order

	order = numpy.asarray(order,dtype=long)
	if order.size == 1:
		order = order * numpy.ones(len(knots),dtype=long)
	
	# the user can pass an arbitrary linear combination of penalty orders
	penorder = [dict() for i in xrange(len(knots))]
	for o,coefficients in penalties.items():
		if int(o) <= 0:
			raise ValueError, "Penalty order must by > 0 (not %s)" % o
		if coefficients is None:
			# if no coefficient is specified, use the smoothness (old behavior)
			coefficients = smooth
		coefficients = numpy.asarray(coefficients,dtype=float)	
		if coefficients.size == 1:
			coefficients = coefficients * numpy.ones(len(knots),dtype=float)
		for i,coeff in enumerate(coefficients):
			penorder[i][int(o)] = coeff
	
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

	def calcP(nsplines, knots, dim, order, porders):
		nspl = nsplines[dim]
		knots = knots[dim]

		def divdiff(knots, m, i):
			# Calculate divided difference coefficients
			# in order to estimate derivatives.

			if m == 0:
				return numpy.asarray([1.])

			num = numpy.append(0,divdiff(knots,m-1,i+1)) - numpy.append(divdiff(knots,m-1,i),0)
			dem = (knots[i+order+1] - knots[i+m])/(order-(m-1))
			return num/dem

		def penalty_matrix(penorder):
			D = numpy.zeros((nspl-penorder,nspl),dtype=float)
			for i in range(0, len(D)):
				D[i][i:i+penorder+1] = divdiff(knots,penorder,i)
			return numpy.asmatrix(D)

		D = penalty_matrix(porders.keys()[0])
		DtD = porders.values()[0] * D.transpose() * D
		
		for porder,coeff in porders.items()[1:]:
			D1 = penalty_matrix(porder)
			DtD += coeff * D1.transpose() * D1

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
	# P = smooth*P

	sidelen = numpy.product(nsplines)
	a = numpy.reshape(numpy.zeros(sidelen,float),nsplines)

	print "Reticulating splines..."

	n = 0
	while n < iterations:
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

		if n > 1:
			# fit for residual on further iterations
			x = numpy.reshape(a,r.shape)
			r = numpy.asarray(r - F*x)
			resid = (r**2).sum()
			print 'The sum of squared residuals is %e'%resid

		result = numpy.linalg.lstsq(F, r)
		
		coefficients = numpy.reshape(result[0],nsplines)
		
		if n == 1:
			a = coefficients
		else:
			a = a + coefficients
			
		
		

	table.coefficients = a
	return table

def monotonize(table,monodim=0):
	"""Use the t-spline hammer to enforce monotonicity along one axis"""
	print "Futzing with t-spline basis"
	
	nsplines = []
	for i in range(0,len(table.knots)):
		if table.periods[i] == 0:
			nsplines.append(len(table.knots[i])-table.order[i]-1)
		else:
			nsplines.append(len(table.knots[i]))
	
	# multiplying by a lower-triangular matrix sums b-spline coefficients
	# to yield t-spline (cumulative) coefficients
	L = numpy.tril(numpy.ones((nsplines[monodim],nsplines[monodim])))
	# the b-spline coefficients are the differences between the t-spline coefficients
	Linv = numpy.linalg.inv(L)
	
	def futz(b_coefficients):
		# first, convert b-spline coefficients to t-spline coefficients
		coefficients = numpy.dot(Linv,b_coefficients)
		for i in xrange(len(coefficients)):
			a = coefficients[i]
			if a < 0:
				print 't-spline coeff %d = %e' % (i,a)
				if (i > 0) and (coefficients[i-1]+a >= 0): 
					# we're coming out of an (over-)ring; add back'erds
					coefficients[i-1] += a
					coefficients[i] = 0
				elif i+1 < len(coefficients):
					# we're going in to an (under-)ring; add forwards
					coefficients[i+1] += a
					coefficients[i] = 0

		# now, convert back to a b-spline basis
		coefficients = numpy.dot(L,coefficients)
		return coefficients
	table.coefficients = numpy.apply_along_axis(futz,monodim,table.coefficients)
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

