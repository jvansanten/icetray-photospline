from . import splinetable
import pyfits
import numpy

def write(table,path):
	"""
	Write a SplineTable to disk as a FITS file.
	
	:param table: the SplineTable to be written
	:param path: where to save the FITS file.
	
	.. warning:: pyfits will fail to write the file if it already exists.
	"""
	data = pyfits.PrimaryHDU(table.coefficients)
	data.header.update('TYPE','Spline Coefficient Table')

	if getattr(table.order,'__iter__',False):
		for i in range(0,len(table.order)):
			data.header.update('ORDER%d' % i,table.order[i],
			    'B-Spline Order')
	else:
		data.header.update('ORDER',table.order,'B-Spline Order')

	for i in range(0,len(table.periods)):
		data.header.update('PERIOD%d' % i,table.periods[i])

	data.header.update('BIAS',table.bias)
	data.header.update('GEOMETRY',table.geometry)
	data.header.update('LEVEL',table.level)
	data.header.update('GEOTYPE',table.geotype)
	data.header.update('NGROUP',table.ngroup)
	data.header.update('PARITY',table.parity)

	hdulist = pyfits.HDUList([data])

	for i in range(0,len(table.knots)):
		knothdu = pyfits.ImageHDU(table.knots[i],name='KNOTS%d' % i)
		hdulist.append(knothdu)
	
	extents = []
	for i in range(0,len(table.knots)):
		if getattr(table.order,'__iter__',False):
			order = table.order[i]
		else:
			order = table.order

		if i < len(table.extents):
			ext = list(table.extents[i])
		else:
			ext = [table.knots[i][order], table.knots[i][-order-1]]
		extents += ext
	extenthdu = pyfits.ImageHDU(numpy.array(extents, dtype=float), name='EXTENTS')
	hdulist.append(extenthdu)

	hdulist.writeto(path)

def read(path, memmap=False):
	"""
	Read a SplineTable from a FITS file on disk
	
	:param path: the filesystem path to read from
	:param memmap: memmap the underlying fits file this causes the underlying file handle to remain open indefinitely
	:returns: SplineTable - the spline surface stored in the given file
	"""
	
	file = pyfits.open(path, memmap=memmap)
	table = splinetable.SplineTable()

	data = file[0]
	table.coefficients = data.data
	table.periods = []
	table.knots = []
	for i in range(0,table.coefficients.ndim):
		try:
			table.periods.append(data.header['PERIOD%d' % i])
		except KeyError:
			table.periods.append(0.)
		table.knots.append(file['KNOTS%d' % i].data)


	try:
		table.order = data.header['ORDER']
		order = [table.order]*table.coefficients.ndim
	except KeyError:
		table.order = []
		for i in range(0,table.coefficients.ndim):
			table.order.append(data.header['ORDER%d' % i])
		order = table.order

	try:
		extents = file['EXTENTS'].data
	except KeyError:
		extents = []
	
	if len(extents) != 2*table.coefficients.ndim:
		extents = []
		for i in range(0,table.coefficients.ndim):
			extents += [table.knots[i][order[i]], table.knots[i][-order[i]-1]]

	table.extents = list(zip(extents[:-1:2], extents[1::2]))

	try:
		table.bias = data.header['BIAS']
	except KeyError:
		pass
	try:
		table.geometry = data.header['GEOMETRY']
	except KeyError:
		pass
	try:
		table.level = data.header['LEVEL']
	except KeyError:
		pass
	try:
		table.geotype = data.header['GEOTYPE']
	except KeyError:
		pass
	try:
		table.ngroup = data.header['NGROUP']
	except KeyError:
		pass
	try:
		table.parity = data.header['PARITY']
	except KeyError:
		pass

	return table

