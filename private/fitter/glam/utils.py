import copy
import numpy
from pyphotonics.photonics_table import photonics_table
import glam,bspline,splinefitstable

class TableSlice(object):
	"""A slice of a photonics table, with spline CDF and PDF evaluates."""
	table_pdf = None
	table_cdf = None
	spline_pdf = None
	spline_cdf = None
	edges = None
	centers = None
	def __init__(self, table, spline, slices, density = 1):
		self.table = table
		self.spline = spline
		self.slices = slices[:table.values.ndim]
		self.density = density
		self.make_grid()
		if len(spline.knots) == 4:
			norm = True
			is_log = False
		else:
			norm = False
			is_log = True
		self.slice(is_log, norm)
		self.eval(is_log)

	def collapse(self):
		"""Collapse non-extended dimensions, returning 1 or 2-d arrays
		   for use with matplotlib.
		"""
		new = copy.copy(self)
		ext_dims = [i for i in xrange(len(self.centers)) if len(self.centers[i]) > 1]
		new.edges = [self.edges[i] for i in ext_dims]
		new.centers = [self.centers[i] for i in ext_dims]
		return new
		
	def flatten(self):
		"""Flatten grid to columns for use with Gnuplot"""
		coords = numpy.meshgrid_nd(*self.centers,**{'lex_order':True})
		coords = numpy.column_stack(map(lambda arr: arr.reshape(arr.size),coords))
		bigdat = numpy.column_stack((coords,self.table_cdf.flatten()))
		if self.table_pdf is not None:
			bigdat = numpy.column_stack((bigdat,self.table_pdf.flatten()))
    		bigdat = numpy.column_stack((bigdat,self.spline_cdf.flatten()))
		if self.spline_pdf is not None:
    			bigdat = numpy.column_stack((bigdat,self.spline_pdf.flatten()))
		return bigdat

	def make_grid(self, density = 1):
		slices = self.slices
		centers = [bins[_slice] for bins,_slice in zip(self.table.bin_centers,slices)]
		for i,sl in enumerate(slices):
			if isinstance(sl,int):
				centers[i] = numpy.array([centers[i]])
			elif self.density > 1: # add extra grid points in between the bin centers
				gaps = numpy.diff(centers[i])
				extras = []
				for j in xrange(1,self.density):
					scale = j/float(self.density)
					extras.append(centers[i][:-1]+scale*gaps)
				centers[i] = numpy.concatenate(tuple([centers[i]] + extras))
				centers[i].sort()
		self.centers = centers

		widths = [bins[_slice] for bins,_slice in zip(self.table.bin_widths,slices)]
		for i,sl in enumerate(slices):
			if isinstance(sl,int):
				widths[i] = numpy.array([widths[i]])
			elif self.density > 1:
				# subdividing the widths is much easier!
				rep = self.density*numpy.ones(widths[i].size, dtype=int)
				rep[-1] = 1
				widths[i] = (widths[i]/self.density).repeat(rep)
		edges = [c - w/2.0 for c,w in zip(centers,widths)]
		edges = [numpy.append(e, c[-1]+w[-1]/2.0) for e,c,w in zip(edges,centers,widths)]
		self.edges = edges

		if len(self.spline.knots) == 4:
			# XXX: evaluate CDF at right edge of the time bin
			self.centers[3] = self.edges[3][1:]

	def slice(self, is_log = False, norm = True):
		timing = len(self.spline.knots) == 4
		table_pdf = None
		table_cdf = None
		slices = self.slices
		table = self.table
		if self.table.normed and self.slices[-1] == slice(None): # already normed
			table_cdf = self.table.values[slices].copy()
			if timing:
				print table_cdf.shape, table.bin_widths[-1].size
				table_pdf = numpy.append([0],
				    numpy.diff(table_cdf, axis=-1))/table.bin_widths[-1]
		elif self.table.normed:
			# already normed, but we're slicing perpendicular to time.
			# skip the PDF.
			table_cdf = self.table.values[slices].copy()
		elif len(self.spline.knots) == 3:
			# abs spline, just sum amplitudes
			tslices = list(self.slices)[:3]
			tslices += [slice(None)]
			bigslice = self.table.values[tslices]
			table_cdf = numpy.sum(bigslice, axis=-1)
		elif timing and self.slices[-1] == slice(None): # we're slicing in time, compute CDF for slice
			table_slice = self.table.values[slices]
			#table_cdf = numpy.cumsum(table_slice*table.bin_widths[3], axis=-1)
			table_cdf = numpy.cumsum(table_slice, axis=-1)
			#table_pdf = table_slice
			table_pdf = table_slice/table.bin_widths[3]
			if norm:
				normslices = [slice(None)]*len(table_cdf.shape)
				normslices[-1] = -1
				newshape = list(table_cdf.shape)
				newshape[-1] = 1
				normval = table_cdf[tuple(normslices)].reshape(tuple(newshape))
				table_cdf /= normval
				table_pdf /= normval
		else:
			# we're slicing perpendicular to time, so slice in the dim-t plane
			# , compute the CDF, then slice in t
			tslices = list(self.slices)
			if timing:
				tslices[-1] = slice(None)
			else:
				tslices += [slice(None)]
			print tslices
			bigslice = self.table.values[tslices]
			table_cdf = numpy.cumsum(bigslice, axis=-1)
			# remove integer indexes, since those dimensions are already gone
			tslices = [t for t in tslices if not isinstance(t,int)]
			# now slice at the time we want
			tslices[-1] = slices[-1]
			nslice = [slice(None)]*table_cdf.ndim
			nslice[-1] = -1
			normval = table_cdf[nslice]
			table_cdf = table_cdf[tslices]
			if timing:
				# careful! table_pdf is a view into table.values
				table_pdf = table.values[slices].copy()
			else:
				table_cdf = normval
			if norm:
				table_cdf /= normval
				table_pdf /= normval

		if self.density > 1: # insert zeros into table values
			expanded_shape = tuple([s + (self.density-1)*(s-1) for s in table_cdf.shape])
			insert_slices = [slice(None,None,self.density)]*len(table_cdf.shape)
			t_cdf = numpy.zeros(expanded_shape)
			t_cdf[insert_slices] = table_cdf
			table_cdf = t_cdf
			if timing:
				t_pdf = numpy.zeros(expanded_shape)
				t_pdf[insert_slices] = table_pdf
				table_pdf = t_pdf
		self.table_cdf = table_cdf
		self.table_pdf = table_pdf
		
	def eval_cdf(self, is_log):
		vals = glam.grideval(self.spline, self.centers)
		shape = vals.shape
		vals = vals[numpy.isfinite(vals)] + self.spline.bias
		shape = tuple([shape[i] for i in xrange(len(shape)) if shape[i] > 1])
		vals = vals.reshape(shape)
		if is_log:
			vals = numpy.exp(vals)
		self.spline_cdf = vals
	def eval_pdf(self, is_log):
		
		splfuncs = [bspline.bspline]*4
		splfuncs[-1] = bspline.bspline_deriv
		deriv_basis = [bspline.splinebasis(self.spline.knots[i], self.spline.order[i],self.centers[i],
                    self.spline.periods[i],spline=splfuncs[i]) for i in range(0,len(self.spline.knots))]
		pdf_vals = glam.grideval(self.spline, self.centers, bases = deriv_basis)
		shape = pdf_vals.shape
		pdf_vals = pdf_vals[numpy.isfinite(pdf_vals)] #+ spline.bias
		if is_log:
			# chain rule!
			pdf_vals *= self.spline_cdf
		shape = tuple([shape[i] for i in xrange(len(shape)) if shape[i] > 1])
		pdf_vals = pdf_vals.reshape(shape)
		self.spline_pdf = pdf_vals
		
		
	def eval(self, is_log = False):
		"""is_log => fit is in log-space"""
		self.eval_cdf(is_log)
		# only calculate PDFs for timing fits
		if (len(self.spline.knots)) == 4:
			self.eval_pdf(is_log)

import itertools as it

def meshgrid_nd(*arrays,**kwargs):
	"""
	Creates coordinate tensors from an arbitrary number of one-dimensional
	coordinate vectors. With 2 arguments, the behavior is identical to
	numpy.meshgrid.
	
	This can be useful, e.g. for evaluating a function that takes ndarrays
	on an N-dimensional without resorting to numpy.vectorize or explicit
	loops.
	
	The dimensions of the returned array are transposed: the coordinate given
	by the first argument varies most quickly, followed by the second, etc.
	This matches the behavior of numpy.meshgrid. To get an array with
	dimensions in lexographical order, pass lex_order=True:
	
	>>> x,y,z=numpy.arange(0,3),numpy.arange(4,6),numpy.arange(7,10)
	>>> X,Y,Z=numpy.meshgrid_nd(x,y,z,lex_order=False)
	>>> X
	array([[[0, 1, 2],
	        [0, 1, 2]],

	       [[0, 1, 2],
	        [0, 1, 2]],

	       [[0, 1, 2],
	        [0, 1, 2]]])
	>>> Y
	array([[[4, 4, 4],
	        [5, 5, 5]],

	       [[4, 4, 4],
	        [5, 5, 5]],

	       [[4, 4, 4],
	        [5, 5, 5]]])
	>>> Z
	array([[[7, 7, 7],
	        [7, 7, 7]],

	       [[8, 8, 8],
	        [8, 8, 8]],

	       [[9, 9, 9],
	        [9, 9, 9]]])
	>>> X,Y,Z=numpy.meshgrid_nd(x,y,z,lex_order=True)
	>>> X
	array([[[0, 0, 0],
	        [0, 0, 0]],

	       [[1, 1, 1],
	        [1, 1, 1]],

	       [[2, 2, 2],
	        [2, 2, 2]]])
	>>> Y
	array([[[4, 4, 4],
	        [5, 5, 5]],

	       [[4, 4, 4],
	        [5, 5, 5]],

	       [[4, 4, 4],
	        [5, 5, 5]]])
	>>> Z
	array([[[7, 8, 9],
	        [7, 8, 9]],

	       [[7, 8, 9],
	        [7, 8, 9]],

	       [[7, 8, 9],
	        [7, 8, 9]]])
	
	"""
	asarrays = map(n.asarray,arrays)
	for ar in asarrays:
		if len(ar.shape) != 1: 
			raise ValueError, "arguments must be 1-d arrays"
	dims = map(len,asarrays)
	out = []
	nD = len(dims)
	for i,arr,dim in it.izip(it.count(),asarrays,dims):
		shape = [1]*nD
		shape[nD-1-i] = dim
		x = arr.reshape(*shape)
		for j,k in it.izip(xrange(nD),reversed(xrange(nD))):
			if k==i: continue
			x = x.repeat(dims[k],axis=j)
		if kwargs.get('lex_order',False): x = x.transpose()
		out.append(x)
	return tuple(out)
	
# hook it in to numpy
numpy.meshgrid_nd = meshgrid_nd

import numpy.core.numeric as _nx
from numpy.core.numeric import asarray, zeros, newaxis, outer, \
	 concatenate, isscalar, array, asanyarray
from numpy.core.fromnumeric import product, reshape
def apply_along_axes(func1d,axis,arrs,*args):
	"""
	Apply a function to 1-D slices along the given axis.

	Execute `func1d(a, *args)` where `func1d` operates on a set of 1-D arrays and `a`
	is a 1-D slice of `arr` along `axis`.

	Parameters
	----------
	func1d : function
		This function should accept 1-D arrays. It is applied to 1-D
		slices of `arr` along the specified axis.
	axis : integer
		Axis along which `arr` is sliced.
	arrs : tuple 
		tuple of input arrays. All arrays must have the same shape
	args : any
		Additional arguments to `func1d`.

	Returns
	-------
	outarr : ndarray
		The output array. The shape of `outarr` is identical to the shape of
		`arr`, except along the `axis` dimension, where the length of `outarr`
		is equal to the size of the return value of `func1d`.  If `func1d`
		returns a scalar `outarr` will have one fewer dimensions than `arr`.

	See Also
	--------
	apply_over_axis : Apply a function over 1-D slices of a single array.
	apply_over_axes : Apply a function repeatedly over multiple axes.

	"""
	arrs = map(asarray,arrs)
	arr = arrs[0]
	nd = arr.ndim
	if axis < 0:
		axis += nd
	if (axis >= nd):
		raise ValueError("axis must be less than arr.ndim; axis=%d, rank=%d."
			% (axis,nd))
	ind = [0]*(nd-1)
	i = zeros(nd,'O')
	indlist = range(nd)
	indlist.remove(axis)
	i[axis] = slice(None,None)
	outshape = asarray(arr.shape).take(indlist)
	for arr in arrs[1:]: 
	if tuple(asarray(arr.shape).take(indlist)) != tuple(outshape):
		raise ValueError("Shape of all input arrays must match in all but the selected dimension.")
	i.put(indlist, ind)
	arglist = tuple(map(lambda arr: arr[tuple(i.tolist())],arrs)) + args
	res = func1d(*arglist)
	#  if res is a number, then we have a smaller output array
	if isscalar(res):
		outarr = zeros(outshape,asarray(res).dtype)
		outarr[tuple(ind)] = res
		Ntot = product(outshape)
		k = 1
		while k < Ntot:
			# increment the index
			ind[-1] += 1
			n = -1
			while (ind[n] >= outshape[n]) and (n > (1-nd)):
				ind[n-1] += 1
				ind[n] = 0
				n -= 1
			i.put(indlist,ind)
			arglist = tuple(map(lambda arr: arr[tuple(i.tolist())],arrs)) + args
			res = func1d(*arglist)
			outarr[tuple(ind)] = res
			k += 1
		return outarr
	else:
		Ntot = product(outshape)
		holdshape = outshape
		outshape = list(arr.shape)
		outshape[axis] = len(res)
		outarr = zeros(outshape,asarray(res).dtype)
		outarr[tuple(i.tolist())] = res
		k = 1
		while k < Ntot:
			# increment the index
			ind[-1] += 1
			n = -1
			while (ind[n] >= holdshape[n]) and (n > (1-nd)):
				ind[n-1] += 1
				ind[n] = 0
				n -= 1
			i.put(indlist, ind)
			arglist = tuple(map(lambda arr: arr[tuple(i.tolist())],arrs)) + args
			res = func1d(*arglist)
			#res = func1d(arr[tuple(i.tolist())],*args)
			outarr[tuple(i.tolist())] = res
			k += 1
		return outarr

numpy.apply_along_axes = apply_along_axes	
