import splinetable
import pyfits

def write(table,path):
	data = pyfits.PrimaryHDU(table.coefficients)
	data.header.update('TYPE','Spline Coefficient Table')
	data.header.update('ORDER',table.order,'B-Spline Order')

	for i in range(0,len(table.periods)):
		data.header.update('PERIOD%d' % i,table.periods[i])

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
	for i in range(0,table.coefficients.ndim):
		table.periods.append(data.header['PERIOD%d' % i])
		table.knots.append(file['KNOTS%d' % i].data)

	table.order = data.header['ORDER']
	return table

