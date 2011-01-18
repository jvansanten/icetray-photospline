import numpy
from icecube.photospline.glam import glam
from icecube.photospline import splinefitstable, glamhist
import sys

# Hard-coded params

bins = 50
nknots = 15
smooth = 0.1

zerooffset = 1e-9

def link(z):
	return numpy.log(z + zerooffset) - numpy.log(zerooffset)

# Real code

data = numpy.loadtxt(sys.argv[1])
table = glamhist.fithist(data, numpy.ones(data.shape[0]),bins, nknots, smooth, link)
print "Saving coefficients to %s..." % (sys.argv[1]+".pspl.fits")
splinefitstable.write(table,sys.argv[1]+".pspl.fits")
