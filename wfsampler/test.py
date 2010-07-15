import hippo
import numpy
import wfsample
from glam import spglam, splinefitstable

#firstcoords = (14, 25, 300, 10)
firstcoords = (30,30,30,130)
nphotons = 50000

app = hippo.HDApp()
cdf = splinefitstable.read("/tmp/normed_cdf_test_lotsaknotsa.fits")
canvas = app.canvas()
das = []

def makeplots(coords,reg=False):
	samp = hippo.DataArray()
	wfsamples = wfsample.wfsample("/tmp/normed_cdf_test_lotsaknotsa.fits", coords[:3], size=nphotons)
	samp['Spline Waveform'] = wfsamples
	samp['Pandel'] = wfsample.pandelsample(numpy.hypot(coords[0],coords[2]),nphotons)
	print reg
	if reg:
		samp.register('Sampled Photons')

	dapdf = hippo.DataArray()
	times = numpy.linspace(0,7000,1000)
	dapdf['Times'] = times
	pdf = numpy.exp(spglam.grideval(cdf, [[coords[0]], [coords[1]], [coords[2]], times, [coords[3]]]).reshape(len(times)))
	pdf = numpy.append([0.],numpy.diff(pdf)/numpy.diff(times))
	dapdf['PDF'] = pdf*nphotons/2.
	if reg:
		dapdf.register('PDF')

	disp = hippo.Display('Histogram', samp, ('Spline Waveform',))
	#dapdf['PDF'] = pdf * disp.getDataRep().getBinWidth['Y'] * disp.getDataRep().
	disp.addDataRep('Histogram', samp.dataSource(), ('Pandel',))
	disp.addDataRep('Strip Chart', dapdf.dataSource(), ('Times', 'PDF'))
	reps = disp.getDataReps()
	reps[1].setColor('red')
	reps[2].setColor('blue')
	#reps[2].normalizeTo(reps[0])
	canvas.addDisplay(disp)

	das.append((samp, dapdf))

makeplots(firstcoords, True)

