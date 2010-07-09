import splinetable
import pyfits

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

	hdulist.writeto(path)

def read(path):
	file = pyfits.open(path)
	table = splinetable.SplineTable()

	data = file[0]
	table.coefficients = data.data
	table.periods = []
	table.knots = []
	for i in range(0,table.coefficients.ndim):
		table.periods.append(data.header['PERIOD%d' % i])
		table.knots.append(file['KNOTS%d' % i].data)

	try:
		table.order = data.header['ORDER']
	except:
		table.order = []
		for i in range(0,table.coefficients.ndim):
			table.order.append(data.header['ORDER%d' % i])

	try:
		table.bias = data.header['BIAS']
	except:
		pass
	try:
		table.geometry = data.header['GEOMETRY']
	except:
		pass

	return table

