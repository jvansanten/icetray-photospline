from glam import splinefitstable
from optparse import OptionParser

from pyphotonics.photonics_table import *

import sys
import os
import numpy

# Hard-coded params

#nknots =[17, 6, 12, 25] # [r, phi, z, t]  For Nathan/Jakob's binning

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
optparser.add_option("-s", "--smooth", dest="smooth", type="float",
             help="smoothness coefficient", default=10.0)
optparser.add_option("--prob", dest="prob", action="store_true",
             help="Fit only the normalized CDFs", default=False)
optparser.add_option("--abs", dest="abs", action="store_true",
             help="Fit only the total amplitude in each cell", default=False)
(opts, args) = optparser.parse_args()

# by default, do both fits
if not opts.prob and not opts.abs:
	opts.prob = opts.abs = True

def check_exists(outputfile):
    if os.path.exists(outputfile):
        if raw_input("File %s exists. Overwrite? (y/n)" % outputfile) == 'y':
            os.unlink(outputfile)
        else:
            sys.exit()

if len(args) < 2:
    abs_outputfile = args[0]+"abs.pspl.fits"
    prob_outfile = args[0]+"prob.pspl.fits"
else:
    abs_outputfile = args[1]+".abs.fits"
    prob_outputfile = args[1]+".prob.fits"

if opts.prob: check_exists(prob_outputfile)
if opts.abs: check_exists(abs_outputfile)

smooth = opts.smooth

# Real code
import spglam as glam

table = photonics_table(args[0])

table.convert_to_level1()

# check for a sane normalization
if (Efficiency.DIFFERENTIAL & table.header['efficiency']):
	raise ValueError, "This appears to be a dP/dt table. Don't do that, okay?"
if (not Efficiency.N_PHOTON & table.header['efficiency']):
	raise ValueError, "This table does not appear to be normalized."

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

print "Core knots:", nknots

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

# use quadratic spacing in time
coreknots[3] = numpy.linspace(0,7000**(0.5),nknots[3])**2
endgap = coreknots[3][-1] - coreknots[3][-2]
tknots = numpy.concatenate(([-10,-5,-3,-1], coreknots[3], 7000 + endgap*numpy.arange(1,8)))

def spline_spec(ndim):
   if ndim > 3:
       order = [2,2,2,3]        # Quadric splines for t to get smooth derivatives
       penalties = {2:[smooth]*3 + [0], # penalize curvature in rho,z,phi
                    3:[0]*3 + [smooth]} # order 3 in time CDF => order 2 in time PDF
       knots = [rknots, thetaknots, zknots, tknots]
   else:
       order = [2,2,2]    # Quadric splines to get smooth derivatives
       penalties = {2:[smooth]*3}    # Penalize curvature 
       # XXX HACK: add more knots near the cascade
       extras = numpy.logspace(-1,1,5)
       zk = numpy.concatenate((zknots[abs(zknots) > 10], -extras, extras, [0]))
       zk.sort()
       extras = numpy.logspace(-1,1,10)
       rk = numpy.concatenate((rknots[(rknots > 10)|(rknots < 0)], extras))
       rk.sort()
       knots = [rk, thetaknots, zk]
   return order, penalties, knots


# Take cumulative sum to get the CDF, and adjust fit points to be
# the right edges of the time bins, where the CDF is measured.
table.values = numpy.cumsum(table.values, axis=3)
table.bin_centers[3] += table.bin_widths[3]/2.

# HACK: Move first and last angular bins to 0 and 180
table.bin_centers[1][0] = 0
table.bin_centers[1][table.bin_centers[1].size - 1] = 180

print "Loaded histogram with dimensions ", table.shape()

norm = table.values[:,:,:,-1]

if opts.abs:
	z = numpy.log(norm)

	# remove NaNs and infinites from the fit
	w = numpy.ones(norm.shape)
	#w = numpy.sum(table.weights, axis=-1)
	w[numpy.logical_not(numpy.isfinite(z))] = 0
	z[numpy.logical_not(numpy.isfinite(z))] = 0

	# XXX HACK: de-weight the first radial bin everywhere
	#w[:3,:,:] = 0 

	order, penalties, knots = spline_spec(3)

	print 'Number of knots used: ',[len(a) for a in knots]
	print "Beginning spline fit for abs table..."
	spline = glam.fit(z,w,table.bin_centers[:3],knots,order,smooth,penalties=penalties)

	print "Saving table to %s..." % abs_outputfile
	splinefitstable.write(spline, abs_outputfile)

if opts.prob:
	z = table.values / norm.reshape(norm.shape + (1,))
	# XXX HACK: ignore weights for normalized timing
	w = numpy.ones(table.weights.shape)
	order, penalties, knots = spline_spec(4)

	print 'Number of knots used: ',[len(a) for a in knots]
	print "Beginning spline fit for timing table..."
	spline = glam.fit(z,table.weights,table.bin_centers,knots,order,smooth,penalties=penalties,monodim=3)

	print "Saving table to %s..." % prob_outputfile
	splinefitstable.write(spline, prob_outputfile)


# smoothed = glam.grideval(spline, table.bin_centers)
# resid = (smoothed - table.values)[table.weights != 0]
# fracresid = ((smoothed - table.values)/table.values)[table.weights != 0]
# 
# 
# print "Fit Statistics:"
# print "\tMaximum Deviation from Data:",numpy.max(numpy.abs(resid))
# print "\tRMS Deviation from Data:",numpy.sqrt(numpy.mean(resid**2))
# print "\tMax Fractional Deviation from Data:",numpy.max(numpy.abs(fracresid))
# print "\tMean Fractional Deviation from Data:",numpy.mean(numpy.abs(fracresid))

