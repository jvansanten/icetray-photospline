__all__=['glam', 'spglam', 'photonics', 'splinetable']

try:
	import numpy
	from icecube.load_pybindings import load_pybindings
	load_pybindings(__name__, __path__)
	del numpy
except ImportError:
	pass

try:
	from . import splinetable
	from . import splinefitstable
	from . import glam
	try:
		import spglam
	except ImportError:
		pass
	from . import numpy_extensions
	from . import photonics
except ImportError:
	pass

