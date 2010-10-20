from setuptools import setup, Feature, Extension
from numpy.distutils.system_info import get_info
import sys, os

import ConfigParser
cp = ConfigParser.ConfigParser()
cp.read('setup.cfg')

if cp.has_option('blas','lapack_libs'):
	lapack_libs = cp.get('blas','lapack_libs').split(',')
else:
	lapack_libs = ['goto2']
	

inc_dirs = []
from numpy.distutils.misc_util import get_numpy_include_dirs
inc_dirs.extend(get_numpy_include_dirs())

photonics_version = '1.67'
if 'I3_PORTS' in os.environ:
	photonics_incdir = os.environ['I3_PORTS'] + (
	    '/include/photonics-%s/' % (photonics_version))
	photonics_libdir = os.environ['I3_PORTS'] + (
	    '/lib/photonics-%s/' % (photonics_version))
	print photonics_libdir
else:
	photonics_incdir = photonics_libdir = "" 

expand = lambda path: os.path.expanduser(os.path.expandvars(path))

if cp.has_option('photonics','path'):
	# We a complete distribution of photonics, including the level2
	# headers. Build with this one.
	photonics_dir = expand(cp.get('photonics','path'))
	photo2numpy = Extension("photo2numpy",
	    sources=["photo2numpy/photo2numpy.c"],
	    include_dirs=inc_dirs + [photonics_dir+'/lib',
	        photonics_dir+'/level2'],
	    extra_objects=[photonics_dir+'/lib/.libs/libphotonics.a',
	        photonics_dir+'/level2/.libs/liblevel2amasim.a'],
	    libraries = ['m'] )
elif os.path.isdir(photonics_libdir):
	# We only have the I3 ports install of photonics, which is
	# missing level2_reader.h. Skip the bits that read l2 tables.
	photonics_libs = ['photonics']
	objs = [photonics_libdir + ('/lib%s.a' % lib) 
	    for lib in photonics_libs]
	photo2numpy = Extension("photo2numpy",
	    sources=["photo2numpy/photo2numpy.c"],
	    include_dirs = inc_dirs + [photonics_incdir],
	    extra_objects = objs,
	    define_macros = [('SKIP_LEVEL2',None)],
	    libraries = ['m'] )
else:
	print "Couldn't find photonics! Skipping build of photo2numpy."
	photo2numpy = None
		
spglam = Extension("spglam",
    sources=["cfitter/glam.c", "cfitter/nnls.c", 
        "cfitter/cholesky_solve.c", "cfitter/splineutil.c", 
        "cfitter/pyglam.c","lib/bspline.c"],
    include_dirs=inc_dirs + ['../public'],
    libraries = ['m', 'cholmod', 'ccolamd', 'colamd', 'amd', 'camd', 
        'metis', 'spqr','gfortran','stdc++'] + lapack_libs,
    undef_macros = ['NDEBUG'] )

spglam_feature = Feature('C version of GLAM (requires SuiteSparse)',
    standard=True, ext_modules=[spglam])

features = {'spglam': spglam_feature}

if photo2numpy is not None:
	photo2numpy_feature = Feature('Photonics table support',
	    standard=True, ext_modules=[photo2numpy])
	features['photo2numpy'] = photo2numpy_feature

def setup_package():
	setup(name = 'photospline', version = '1e-12',
	    author = 'Nathan Whitehorn, Jakob van Santen, Sven Lafebre',
	    package_dir = {'glam': 'fitter/glam',
	         'pyphotonics': 'fitter/pyphotonics'},
	    packages = ["glam","pyphotonics"],
	    features = features
	)

if __name__ == "__main__":
	setup_package()
