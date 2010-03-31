from glam import splinefitstable
from optparse import OptionParser

from pyphotonics.photonics_table import *

import sys
import os
import numpy

# Hard-coded params

#nknots =[17, 6, 12, 25] # [r, phi, z, t]  For Nathan/Jakob's binning
smooth = 10.

# Parse arguments

usage = "usage: %prog [options] table.pt [output.fits]"
optparser = OptionParser(usage=usage)
optparser.add_option("-e", "--epsilon", dest="epsilon", type="float",
             help="epsilon with which to replace table zeros")
optparser.add_option("-r", "--rknots", dest="rknots", type="int",
             help="number of knots in radial dimension")
optparser.add_option("-f", "--fknots", dest="fknots", type="int",
             help="number of knots in angular dimension")
optparser.add_option("-z", "--zknots", dest="zknots", type="int",
             help="number of knots in longitudinal dimension")
optparser.add_option("-t", "--tknots", dest="tknots", type="int",
             help="number of knots in time dimension")
(opts, args) = optparser.parse_args()

# Real code
from glam import spglam as glam

table = photonics_table(args[0])

table.convert_to_level1()

nknots = [17, 6, 12]
if table.ndim() > 3:
    nknots.append(25) # [t]

if opts.rknots:
    nknots[0] = opts.rknots
if opts.fknots:
    nknots[1] = opts.fknots
if opts.zknots:
    nknots[2] = opts.zknots
if opts.tknots and table.ndim() > 3:
    nknots[3] = opts.tknots

print "Knots used to fit table:", nknots

# Compute knot locations using a sparsified version of the bin centers as
#    central knots.

radial_extent = 600

coreknots = [table.bin_centers[i][numpy.unique(numpy.int32(numpy.linspace(0,len(table.bin_centers[i])-1,nknots[i])))] for i in range(0,table.ndim())]

# Now append the extra knots off both ends of the axis in order to provide
# full support at the boundaries

rknots     = numpy.append(numpy.append([-1, -0.5, -0.1], coreknots[0]),
                          numpy.asarray([100, 200, 300, 400, 500]) + radial_extent)
thetaknots = numpy.append(numpy.append([-1, -0.5, -0.1], coreknots[1]),
                          [180.1,180.2,180.3,180.4,180.5])
zknots     = numpy.append(numpy.append([-800, -700, -600], coreknots[2]),
                          [600,700,800,900,1000])

# Use log spacing for time
if table.ndim() > 3:
    coreknots[3] = numpy.logspace(0,numpy.log10(7000),nknots[3])
    tknots = numpy.append(numpy.append([-1,-0.5,-0.25,0], coreknots[3]),
                          [7100, 7150, 7200, 7300, 7400])

if table.ndim() > 3:
    order = [2,2,2,3]        # Quadric splines for t to get smooth derivatives
    penorder = [2,2,2,1]    # Penalize non-constancy in the CDF
    knots = [rknots, thetaknots, zknots, tknots]
else:
    order = [2,2,2]    # Quadric splines for t to get smooth derivatives
    penorder = [2,2,2]    # Penalize non-constancy in the CDF
    knots = [rknots, thetaknots, zknots]

print 'Number of knots used: ',[len(a) for a in knots]

if opts.epsilon != None:
    table.weights[table.values == 0] = 0.01
    table.values = table.values + opts.epsilon

# Take cumulative sum to get the CDF, and adjust fit points to be
# the right edges of the time bins, where the CDF is measured.
#if table.ndim() > 3:
    #table.values = numpy.cumsum(table.values, axis=3)
    #table.bin_centers[3] += table.bin_widths[3]/2.

#if len(table.bin_centers) > table.ndim():
    #table.bin_centers = table.bin_centers[0:table.ndim()]

# HACK: Move first and last angular bins to 0 and 180
table.bin_centers[1][0] = 0
table.bin_centers[1][table.bin_centers[1].size - 1] = 180

# Convert the input to log-space and drop any NaNs or infinites from the fit
table.values = numpy.log(table.values)

table.remove_nans_and_infinites()

print "Loaded histogram with dimensions ", table.shape()

if len(args) < 2:
    outputfile = args[0]+".pspl.fits"
else:
    outputfile = args[1]

if os.path.exists(outputfile):
    if raw_input("File %s exists. Overwrite? (y/n)" % outputfile) == 'y':
        os.unlink(outputfile)
    else:
        sys.exit()

print "Beginning spline fit..."
spline = glam.fit(table.values,table.weights,table.bin_centers,knots,order,smooth,penorder=penorder)

print "Saving table to %s..." % outputfile
splinefitstable.write(spline, outputfile)

smoothed = glam.grideval(spline, table.bin_centers)
resid = (smoothed - table.values)[table.weights != 0]
fracresid = ((smoothed - table.values)/table.values)[table.weights != 0]

print "Fit Statistics:"
print "\tMaximum Deviation from Data:",numpy.max(numpy.abs(resid))
print "\tRMS Deviation from Data:",numpy.sqrt(numpy.mean(resid**2))
print "\tMax Fractional Deviation from Data:",numpy.max(numpy.abs(fracresid))
print "\tMean Fractional Deviation from Data:",numpy.mean(numpy.abs(fracresid))
