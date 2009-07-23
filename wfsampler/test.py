import hippo
import numpy
import wfsample
from glam import spglam, splinefitstable

coords = (14, 25, 300, 10)
nphotons = 5000

app = hippo.HDApp()

samp = hippo.DataArray()
samp['Spline Waveform'] = wfsample.wfsample("dndt.fits", coords, size=nphotons)
samp['Pandel'] = wfsample.pandelsample(300,nphotons)
samp.register('Sampled Photons')

dndt = splinefitstable.read("dndt.fits")
dapdf = hippo.DataArray()
times = numpy.linspace(0,7000,1000)
dapdf['Times'] = times
dapdf['PDF'] = numpy.exp(spglam.grideval(dndt, [[coords[0]], [coords[1]], [coords[2]], times, [coords[3]]]).reshape(len(times)))*1e16
dapdf.register('PDF')

canvas = app.canvas()
disp = hippo.Display('Histogram', samp, ('Spline Waveform',))
disp.addDataRep('Histogram', samp.dataSource(), ('Pandel',))
disp.addDataRep('Strip Chart', dapdf.dataSource(), ('Times', 'PDF'))
reps = disp.getDataReps()
reps[1].setColor('red')
reps[2].setColor('blue')
canvas.addDisplay(disp)
