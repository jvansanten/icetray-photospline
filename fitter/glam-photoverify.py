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

usage = "usage: %prog [options] table.pt table.fits"
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
    tree_d = True

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

spline   = splinefitstable.read(args[1])
smoothed = glam.grideval(spline,table.bin_centers)

print "Loaded spline"

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
bigdat = numpy.column_stack((bigdat,smoothed.reshape(smoothed.size)))
bigdat = bigdat[numpy.isfinite(bigdat[:,table.ndim()+1])]

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
        if three_d:
            raw  = Gnuplot.Data(sample[:,xdim], sample[:,ydim], sample[:,table.ndim()],   title="Raw")
            fit = Gnuplot.Data(sample[:,xdim], sample[:,ydim], sample[:,table.ndim()+1], title="Fit")
            gp.ylabel(axis_labels[ydim])
            gp.zlabel("log Photon density")
            gp.splot(raw,fit)
        else:
            raw = Gnuplot.Data(sample[:,xdim], sample[:,table.ndim()],   title="Raw")
            fit = Gnuplot.Data(sample[:,xdim], sample[:,table.ndim()+1], title="Fit")
            gp.ylabel("log Photon density")
            gp.plot(raw,fit)

        printdiffstats(sample[:,table.ndim()],sample[:,table.ndim()+1])
        raw_input('Press Enter to continue')
    except AssertionError:
        print 'Skipping   slice at  (no raw data)' % i
