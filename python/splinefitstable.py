import splinetable
import pyfits
import numpy

def write(table,path):
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

def read(path):
	file = pyfits.open(path)
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

	table.extents = zip(extents[:-1:2], extents[1::2])

	try:
		table.bias = data.header['BIAS']
	except KeyError:
		pass
	try:
		table.geometry = data.header['GEOMETRY']
	except KeyError:
		pass

	return table

