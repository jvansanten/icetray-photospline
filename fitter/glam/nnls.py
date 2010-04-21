import numpy as n
import scipy.linalg

def lusolve(A,b):
	lu = scipy.linalg.lu_factor(A)
	solution = scipy.linalg.lu_solve(lu,b)
	return solution
	
def cholsolve(A,b):
	chol = scipy.linalg.cho_factor(A)
	solution = scipy.linalg.cho_solve(chol,b)
	return solution


def nnls(A,b,verbose=True):
	"""A mockup of the Lawson/Hanson active set algorithm
	
	See:
	A Comparison of Block Pivoting and Interior-Point Algorithms for Linear Least Squares Problems with Nonnegative Variables
	Author(s): Luis F. Portugal, Joaquim J. Judice, Luis N. Vicente
	Source: Mathematics of Computation, Vol. 63, No. 208 (Oct., 1994), pp. 625-643
	Published by: American Mathematical Society
	Stable URL: http://www.jstor.org/stable/2153286
	"""
	# let's have the proper shape
	A = n.asarray(A)
	b = n.asarray(b).reshape(b.size)
	# step 0
	F = [] # passive (solved for) set
	G = range(A.shape[0]) # active (clamped to zero) set
	x = n.zeros(A.shape[0])
	
	y = -n.dot(A.transpose(),b)
	
	if verbose:
		def log(mesg):
			print mesg
	else:
		def log(mesg):
			pass
	
	iterations = 0
	lstsqs = 0
	while True:
		iterations += 1
		# step 1
		if len(G) == 0:
			log("Active set empty, terminating after %d iterations (%d least squares computed)" % (iterations,lstsqs))
			break # the active set is the whole set, we're done
		r_G = y[G].argmin()
		r = G[r_G]
		# print x,y
		if y[r] >= 0:
			log("Dual vector is all positive, terminating after %d iterations (%d least squares computed)" % (iterations,lstsqs))
			break # x is the optimal solution, we're done
		log("Moving %d into active set" % r)
		F.append(r); F.sort()
		G.remove(r)
		feasible = False
		while not feasible:
			# print 'F:',F
			# print 'G:',G
			# step 2
			log("Unconstrained solve: %s, %s" % (A[:,F].shape,b.shape))
			x_F = n.linalg.lstsq(A[:,F],b)[0]
			lstsqs += 1
			if (x_F >= 0).all():
				x[F] = x_F
				feasible = True
			else:
				# if the new trial solution gained a negative element
				mask = (x_F <= 0)
				theta = x[F]/(x[F] - x_F)
				
				r_F = theta[mask].argmin()
				alpha = theta[mask][r_F]
				r = n.array(F)[mask][r_F]
				x[F] = x[F] + alpha*(x_F-x[F])
				log("Moving %d to passive set" % r)
				F.remove(r)
				G.append(r); G.sort()
		# step 3
		y[:] = 0
		y[G] = n.dot(A[:,G].transpose(),(n.dot(A[:,F],x[F])-b))
	return x

def nnls_normal(AtA,Atb,verbose=True):
	"""A mockup of the Lawson/Hanson active set algorithm for pre-formulated normal equations
	
	This version starts from the unconstrained solution (which may be moderately faster)"""
	# let's have the proper shape
	AtA = n.asarray(AtA)
	Atb = n.asarray(Atb).reshape(Atb.size)
	nvar = AtA.shape[0]
	maxiter = 3*nvar
	# step 0
	F = [] # passive (solved by unconstrained least squares) set
	G = range(nvar) # active (clamped to zero) set
	x = n.zeros(nvar)

	if verbose:
		def log(mesg):
			print mesg
	else:
		def log(mesg):
			pass

	# variant: initialize with unconstrained solution
	x = lusolve(AtA,Atb)
	mask = x < 0
	indices = n.arange(x.size)
	F = list(indices[n.logical_not(mask)])
	G = list(indices[mask])
	x[mask] = 0
	y = n.zeros(x.size)
	y[G] = n.dot(AtA[:,F][G,:],x[F]) - Atb[:,G]
	log(F)
	log(G)
	log(x)
	log(y)

	iterations = 0
	lstsqs = 0
	
	while True:
		iterations += 1
		# step 1
		if iterations > 1:
			if len(G) == 0:
				log("Active set empty, terminating after %d iterations (%d LU solves)" % (iterations,lstsqs+1))
				break # the passive set is the whole set, we're done
			r_G = y[G].argmin()
			r = G[r_G]
			# print x,y
			if y[r] >= 0:
				log("Dual vector is all positive, terminating after %d iterations (%d LU solves)" % (iterations,lstsqs+1))
				break # x is the optimal solution, we're done
			log("Moving %d into active set" % r)
			F.append(r); F.sort()
			G.remove(r)
		feasible = False
		while not feasible:
			# print 'F:',F
			# print 'G:',G
			# step 2
			
			# x_F = n.linalg.lstsq(A[:,F],b)[0]
			# select only the bits of A^T*A that apply to coefficients F
			AtA_F = AtA[:,F][F,:]
			Atb_F = Atb[:,F]
			log("Unconstrained solve: %s, %s" % (AtA_F.shape,Atb_F.shape))
			x_F = lusolve(AtA_F,Atb_F)
			lstsqs += 1
			if (x_F >= 0).all():
				x[F] = x_F
				feasible = True
			else:
				# if the new trial solution gained a negative element,
				# find the worst offending coefficient and move it back to the passive set
				mask = (x_F <= 0)
				theta = x[F]/(x[F] - x_F)
				r_F = theta[mask].argmin()
				alpha = theta[mask][r_F]
				r = n.array(F)[mask][r_F]
				x[F] = x[F] + alpha*(x_F-x[F])
				log("Moving %d to passive set" % r)
				F.remove(r)
				G.append(r); G.sort()
		# step 3
		y[:] = 0
		y[G] = n.dot(AtA[:,F][G,:],x[F]) - Atb[:,G]
	return x

def nnls_normal_block(AtA,Atb,verbose=True):
	"""A mockup of the Portugal/Judice/Vicente block-pivoting algorithm for pre-formulated normal equations
	
	See:
	A Comparison of Block Pivoting and Interior-Point Algorithms for Linear Least Squares Problems with Nonnegative Variables
	Author(s): Luis F. Portugal, Joaquim J. Judice, Luis N. Vicente
	Source: Mathematics of Computation, Vol. 63, No. 208 (Oct., 1994), pp. 625-643
	Published by: American Mathematical Society
	Stable URL: http://www.jstor.org/stable/2153286
	"""
	# let's have the proper shape
	AtA = n.asarray(AtA)
	Atb = n.asarray(Atb).reshape(Atb.size)
	nvar = AtA.shape[0]
	maxiter = 3*nvar
	
	if verbose:
		def log(mesg):
			print mesg
	else:
		def log(mesg):
			pass
	
	# step 0
	F = [] # passive (solved by unconstrained least squares) set
	G = range(nvar) # active (clamped to zero) set
	x = n.zeros(nvar)
	y = -Atb
	
	ninf = nvar + 1 # number of infeasible coefficients
	max_trials = 10 # number of block pivots to try before resorting to Murty's method
	
	p = max_trials
	iterations = 0
	
	while iterations < maxiter:
		iterations += 1
		if (x[F] >= 0).all() and (y[G] >= 0).all():
			log('All coefficients are positive, terminating after %d iterations' % iterations)
			break
		H1 = n.array(F)[x[F] < 0]
		H2 = n.array(G)[y[G] < 0]
		current_ninf = len(H1) + len(H2)
		if current_ninf < ninf:
			ninf = current_ninf
			p = max_trials
		elif current_ninf >= ninf:
			if p >= 1:
				p -= 1
			else: # Murty's method (pick the last infeasible coordinate)
				rmax1 = max(H1)
				rmax2 = max(H2)
				if rmax1 > rmax2:
					H1 = [rmax1]; H2 = []
				else:
					H1 = []; H2 = [rmax2]
		# shuffle infeasible coefficients between sets
		log('infeasibles: %d'%ninf)
		for r in H1:
			F.remove(r); G.append(r)
		for r in H2:
			G.remove(r); F.append(r)
		F.sort(); G.sort()
		AtA_F = AtA[:,F][F,:]
		Atb_F = Atb[:,F]
		log("Unconstrained solve for %d of %d coefficients" % (len(F),nvar))
		x[F] = cholsolve(AtA_F,Atb_F)
		x[G] = 0
		y[F] = 0
		y[G] = n.dot(AtA[:,F][G,:],x[F]) - Atb[:,G]
	if iterations == maxiter:
		print 'Hooo boy, this turned out badly'
	return x
		   
def test(size=5):
	import pylab as p
	# A = n.random.uniform(size=size*size,low=-1,high=1).reshape((size,size))
	A = n.eye(size)
	x_true = n.random.uniform(size=size,low=-1,high=1)
	b = n.dot(A,x_true)
	x = nnls(A,b)
	p.figure()
	p.plot(x_true,drawstyle='steps-post',label='true')
	p.plot(x,drawstyle='steps-post',label='fit')
	p.legend()
	p.show()