from glam import glam, splinetable, splinefitstable
from pyphotonics.photonics_table import *

from optparse import OptionParser

from numpy_extensions import *

import sys
import numpy
import Gnuplot

axis_vars   = ("r", "theta", "z", "t")
axis_labels = ("Perpendicular distance (m)",
               "Source-observer angle (degree)",
               "Parallel distance (m)",
               "Time (ns)")

usage = "usage: %prog [options] table.pt table.fits [table2.fits [ ... ] ]"
optparser = OptionParser(usage=usage)
optparser.add_option("-0", "--dim0", dest="dim0", type="string",
             help="How to plot histogram dimension 0 [x|y|i|<bin>]")
optparser.add_option("-1", "--dim1", dest="dim1", type="string",
             help="How to plot histogram dimension 1 [x|y|i|<bin>]")
optparser.add_option("-2", "--dim2", dest="dim2", type="string",
             help="How to plot histogram dimension 2 [x|y|i|<bin>]")
optparser.add_option("-3", "--dim3", dest="dim3", type="string",
             help='''How to plot histogram dimension 3 [x|y|i|<bin>]                               
                     x    : use this dimension as x variable in plot                               
                     y    : use this dimension as y variable in plot                               
                     i    : iterate over this dimension                                            
                     <bin>: slice at this bin number''')
(opts, args) = optparser.parse_args()

# Load oroginal table
table = photonics_table(args[0])

print "Loaded histogram with dimensions ",table.shape()

three_d = False
xdim = None
ydim = None
idim = None
axes = [opts.dim0, opts.dim1, opts.dim2, opts.dim3]
free_axes = range(table.ndim())

if 'y' in axes:
    ydim = axes.index('y') 
    if ydim >= table.ndim():
        print ydim,"-> y: Table only has dimensions", range(table.ndim())
        sys.exit()
    free_axes.remove(ydim)
    three_d = True

if 'x' in axes:
    xdim = axes.index('x')
    if xdim >= table.ndim():
        print xdim,"-> x: Table only has dimensions", range(table.ndim())
        sys.exit()
else:
    xdim = max(free_axes)
free_axes.remove(xdim)

if 'i' in axes:
    idim = axes.index('i')
    if idim >= table.ndim():
        print idim,"-> i: Table only has dimensions", range(table.ndim())
        sys.exit()
else:
    idim = max(free_axes)
free_axes.remove(idim)

#print 'x:', xdim, '   y:', ydim, '   i:', idim, 'free:', free_axes

# Convert the input to log-space and drop any NaNs or infinites from the fit
table.weights[numpy.logical_not(numpy.isfinite(table.values))] = 0
table.values [numpy.logical_not(numpy.isfinite(table.values))] = 0
table.values = numpy.log(table.values)

# Load the spline fit

smoothed = []
for file in args[1:]:
    spline = splinefitstable.read(file)
    vals   = glam.grideval(spline, table.bin_centers)
    vals   = vals[numpy.isfinite(vals)]
    smoothed.append(vals)
    del spline

print "Loaded", len(smoothed),
if len(smoothed) == 1:
    print "spline file"
else:
    print "spline files"

# ---------------------------------------------------------------------------

def printdiffstats(a,b):
    print "Fit Statistics:"
    print "\tMaximum Deviation from Data:",numpy.max(numpy.abs(a - b))
    print "\tRMS Deviation from Data:",numpy.sqrt(numpy.mean((a - b)**2))
    print "\tMax Fractional Deviation from Data:",numpy.max(numpy.abs((a - b)/b))
    print "\tMean Fractional Deviation from Data:",numpy.mean(numpy.abs((a - b)/b))

print "Building data set..."

coords = n.meshgrid_nd(*table.bin_centers,**{'lex_order':True})
coords = numpy.column_stack(map(lambda arr: arr.reshape(arr.size),coords))
zvec = table.values.reshape(table.values.size)
bigdat = numpy.column_stack((coords,zvec))
for spline in smoothed:
    bigdat = numpy.column_stack((bigdat, spline.reshape(spline.size)))

gp = Gnuplot.Gnuplot()

for i in table.bin_centers[idim]:
    title = "Slice at %s = %f" % (axis_vars[idim], i)
    try:
        sample = bigdat[(bigdat[:,idim] == i)]
        for d in free_axes:
            if axes[d] is not None:
                bin = int(axes[d])
            else:
                bin = int(len(table.bin_centers[d])/2.)
            sample = sample[(sample[:,d] == table.bin_centers[d][bin])] 
            title += ", %s = %f" % (axis_vars[d], table.bin_centers[d][bin])

        gp.xlabel(axis_labels[xdim])
        gp.title(title)
        plots = []
        if three_d:
            pass
            for plot in range(table.ndim()+1,sample.shape[1]):
                plots.append(Gnuplot.Data(sample[:,xdim],
                                          sample[:,ydim],
                                          sample[:,plot],
                                          title="Fit %d" % (plot - table.ndim())))
            plots.append(Gnuplot.Data(sample[:,xdim], sample[:,ydim], sample[:,table.ndim()],   title="Raw"))
            gp.ylabel(axis_labels[ydim])
            gp.zlabel("log Photon density")
            #gp.set_range("zrange", "[-50:0]")
            gp.splot(*plots)
        else:
            for plot in range(table.ndim()+1,sample.shape[1]):
                plots.append(Gnuplot.Data(sample[:,xdim],
                                          sample[:,plot],
                                          title="Fit %d" % (plot - table.ndim())))
            plots.append(Gnuplot.Data(sample[:,xdim], sample[:,table.ndim()],   title="Raw"))
            gp.ylabel("log Photon density")
            gp.plot(*plots)

        printdiffstats(sample[:,table.ndim()],sample[:,table.ndim()+1])
        raw_input('Press Enter to continue')
    except AssertionError:
        print 'Skipping   slice at  (no raw data)' % i
