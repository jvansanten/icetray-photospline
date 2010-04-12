from distutils.core import setup, Extension
import os

if not os.environ.has_key('PHOTONICS_DIR'):
	print "Couldn't find photonics! Set PHOTONICS_DIR to the path to your photonics source."
	exit(1)
else:
	photonics_dir = os.environ['PHOTONICS_DIR'] 

inc_dirs = [photonics_dir+'/lib',photonics_dir+'/level2']
statics = [photonics_dir+'/lib/.libs/libphotonics.a',photonics_dir+'/level2/.libs/liblevel2amasim.a']
libs = ['m']

from numpy.distutils.misc_util import get_numpy_include_dirs
inc_dirs.extend(get_numpy_include_dirs())

setup(
    ext_modules = [
        Extension("photo2numpy", sources=["photo2numpy.c"],
                            include_dirs=inc_dirs,
                            extra_objects=statics,
                            libraries=libs
                            )
    ]
)
