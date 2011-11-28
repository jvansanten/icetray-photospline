__all__=['glam', 'spglam', 'photo2numpy', 'pyphotonics', 'splinetable']

try:
	import numpy
	from icecube.load_pybindings import load_pybindings
	load_pybindings(__name__, __path__)
	del numpy
except ImportError:
	pass

try:
	import splinetable
	import splinefitstable
	import glam
	try:
		import spglam
	except ImportError:
		pass
	import numpy_extensions
	import photo2numpy
	import pyphotonics
except ImportError:
	pass

