import numpy as n
import pylab as p

from glam import bspline

def divdiff(x, y):
	if len(y) == 1:
		return y[0]
	else:
		return (divdiff(x[:-1],y[:-1]) - divdiff(x[1:],y[1:]))/(x[-1] - x[0])

def factorial(x):
	out = 1.0
	num = x
	while num > 1:
		out *= num
		num -= 1
	return out

def make_a_spline(knots, order, z):
	"""Does the divided-differences-of-truncated-power-functions representation make sense?"""
	k = order + 1
	
	vals = n.zeros(z.shape)
	for i, zi in enumerate(z):
		tpf = (zi-knots)
		mask = tpf > 0
		tpf[mask] = tpf[mask]**(k-1)
		tpf[n.logical_not(mask)] = 0
		vals[i] = k*divdiff(knots, tpf)
	return vals

def cbar_simple(x, y, z, bags):
	"""A literal implementation of the local blossom of the convolution"""
	
	fun_x = n.zeros(x.shape)
	for i, xi in enumerate(x):
		fun_y = n.zeros(y.shape)
		mask = xi + y - z > 0
		det = 1
		for bag in bags:
			det *= (xi + y - bag)
		fun_y[mask] = det[mask]
		fun_x[i] = divdiff(y, fun_y)
	return divdiff(x, fun_x)
	
def convolved_coefficient(x, xorder, y, yorder, rho, i):
	k = xorder + 1
	q = yorder + 1
	
	ndeg = k + q - 1
	
	bundle = rho[i:i+ndeg+2]
	xknots = x
	yknots = y
	
	scale = (factorial(k)*factorial(q)/factorial(k+q-1))*(bundle[-1] - bundle[0])/(ndeg + 1)
	bags =  bundle[1:-1]
	cbar = cbar_simple(xknots, yknots, bundle[0], bags)
	if k % 2 != 0:
		cbar *= -1
	
	# integral of a De Boor-normalized spline
	norm = (rho[i+ndeg+1] - rho[i])/(ndeg + 1)
	
	return scale*cbar/norm
	
def blossom(x, y, z, bags):
	pass

def test():
	
	k = 3
	
	x = n.logspace(-1,n.log10(k),k+1)
	
	# x = n.arange(k+1, dtype=float)
	y = 0.8*n.array([-2,-1, 0, 1, 2])[1:-1]
	# y = n.arange(-q/2,q/2 + 1, dtype=float)
	
	rho = n.unique(n.array([xi + yi for xi in x for yi in y]))
	rho.sort()	
	
	k = len(x)-1
	q = len(y)-1
	
	z = n.linspace(min(x.min(),y.min()), max(x.max(),y.max()), 500)
	
	convorder = k + q - 1
	print 'order %d * order %d -> order %d' % (k-1, q-1, convorder)
	# convolution product is a spline of degree k+q
	nsplines = rho.size - convorder - 1
	print nsplines,'splines'
	coeffs = []
	for i in xrange(nsplines):
		coeff = convolved_coefficient(x, k-1, y, q-1, rho, i)
		coeffs.append(coeff)
	coeffs = n.array(coeffs)
	print coeffs
	
	basis = n.asarray(bspline.splinebasis(rho, convorder, z))
	
	evaluate = n.zeros(z.shape)
	
	for i,zi in enumerate(z):
		tpf_x = n.zeros(x.shape)
		for j,xi in enumerate(x):
			# tpf_y = xi + y - zi
			tpf_y = zi - xi + y
			mask = tpf_y > 0
			tpf_y[mask] = tpf_y[mask]**(k+q-1)
			tpf_y[n.logical_not(mask)] = 0
			tpf_x[j] = divdiff(y, tpf_y)
		
		operated = divdiff(x, tpf_x)
		evaluate[i] = (factorial(k)*factorial(q)/factorial(k+q-1))*operated
	
	p.figure()
	p.plot(z, make_a_spline(x, k-1, z), label='f')
	p.plot(z, make_a_spline(y, q-1, z), label='g')
		
	spliff = n.dot(coeffs, basis.transpose())
		
	p.plot(z, spliff, label='f*g, spline expansion')

	p.legend()	