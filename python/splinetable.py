class SplineTable:
	# Context information for the fit
	order = 2
	knots = []
	periods = []

	# The tensor-product basis function coefficients
	coefficients = None

	# logarithmic bias
	bias = 0

	# extent of supported region
	extents = []

	# geometry type
	geometry = 2
	# geometry class of photonics table (0 is point souce, 1 is inf. muon)
	geotype = 0

	# level of photonics table
	level = 1

	# group phase velocity if spline is fit
	# from a photonics table

	ngroup = -1.0
