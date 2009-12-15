import numpy
from glam import splinetable, splinefitstable
from glob import glob
import re, os, sys
import copy

#os.chdir(sys.argv[1])
sourcedir = sys.argv[1]
tables = [(i, re.match(".*_z(-?\d+)_a(\d+).*", i).groups()) for i in glob("%s*.fits"%sourcedir)]

# Convert tables to a useful form and sort it

tables = sorted([(i[0], (int(i[1][0]), int(i[1][1]))) for i in tables], key=lambda tab: tab[1])


# Read in all the actual tables

print 'Table list acquired, reading in tables...'

tables = [(splinefitstable.read(i[0]), i[1]) for i in tables]

print '\tDONE'

def stack_tables(tablist):
	# We expect an array of (splinetable, coordinate) tuples

	bigtab = None

	for table in tablist:
		slice = table[0]
		position = table[1]

		slice.coefficients = slice.coefficients.reshape( \
		    slice.coefficients.shape + (1,))
		knotloc = position + 0.5*(slice.order - 1.0)
		ndim = slice.coefficients.ndim
		if bigtab is None:
			bigtab = slice
			bigtab.knots.append([knotloc])
			bigtab.periods.append(0)
		else:
			bigtab.knots[ndim - 1].append(knotloc)
			bigtab.coefficients = numpy.concatenate(
			    (bigtab.coefficients, slice.coefficients),
			    ndim - 1)

	# Shift the knots (see bsplineinterp.py)
	baseknots = bigtab.knots[ndim - 1]
	baseknots = baseknots + (numpy.max(baseknots)-numpy.min(baseknots))/(2.0*len(baseknots))*(bigtab.order-1)
	interpknots = []
	for i in range (bigtab.order,0,-1):
		interpknots.append(baseknots[0] - i*(baseknots[1] - baseknots[0]))
	interpknots.extend(baseknots)
	interpknots.append(interpknots[len(interpknots)-1] + (interpknots[len(interpknots)-1] - interpknots[len(interpknots)-2]))
	bigtab.knots[ndim - 1] = numpy.asarray(interpknots)

	return bigtab

zpos = numpy.unique([i[1][0] for i in tables])

intermedtables = []

for z in zpos:
	print 'Stacking tables at z =', z
	
	# Select all the tables at this z
	sublist = filter(lambda tab: tab[1][0] == z, tables)
	# Reformat to just one coordinate for stacking
	sublist = [(tab[0], tab[1][1]) for tab in sublist]
	# extend angular range by mirroring next-to-last
	# angle bins (e.g. 170 and 10 deg) to the outside
	# (e.g. 190 and -10) so that 0 and 180 will have
	# support
	print '\textending angular range by hand...'
	lowmirror = [(copy.deepcopy(sublist[1][0]),-sublist[1][1])]
	highmirror = [(copy.deepcopy(sublist[-2][0]),sublist[-1][1]+(sublist[-1][1]-sublist[-2][1]))]
	sublist = lowmirror + sublist + highmirror

	intermedtables.append((stack_tables(sublist), z))

# We no longer need to original tables
del tables

if len(zpos) > 1:
	finaltab = stack_tables(intermedtables)
else:
	finaltab = intermedtables[0][0]

splinefitstable.write(finaltab, '../' + os.path.basename(os.path.abspath(sourcedir)) + '.fits')
