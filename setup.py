from distutils.core import setup, Extension
from numpy.distutils.system_info import get_info
import os

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

def setup_extensions():
	if cp.has_option('photonics','path'):
		photonics_dir = cp.get('photonics','path')
		photo2numpy = Extension("photo2numpy", sources=["photo2numpy/photo2numpy.c"],
							include_dirs  = inc_dirs + [photonics_dir+'/lib',
											 photonics_dir+'/level2'],
							extra_objects = [photonics_dir+'/lib/.libs/libphotonics.a',
											 photonics_dir+'/level2/.libs/liblevel2amasim.a'],
							libraries	  = ['m']
							)
	else:
		print "Couldn't find photonics! Skipping build of photo2numpy..."
		photo2numpy = None
	
	spglam = Extension("spglam", sources = ["cfitter/glam.c","cfitter/splineutil.c","cfitter/pyglam.c","lib/bspline.c"],
						include_dirs = inc_dirs + ['lib'],
						libraries = ['m','cholmod','ccolamd','colamd','amd','spqr','gfortran','stdc++','gfortranbegin'] + lapack_libs,
						)
	
	if photo2numpy is not None:
		return [photo2numpy, spglam]
	else:
		return [spglam]

def setup_package():
	setup(
		name = 'photospline',
		version = '1e-12',
		author = 'Nathan Whitehorn, Jakob van Santen, Sven Lafebre',
		package_dir = {'glam': 'fitter/glam', 'pyphotonics': 'fitter/pyphotonics'},
		packages = ["glam","pyphotonics"],
		
		ext_modules = setup_extensions
	)

if __name__ == "__main__":
	setup_package()