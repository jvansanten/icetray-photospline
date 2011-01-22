class SplineTable:
	# Context information for the fit
	order = 2
	knots = []
	periods = []

	# The tensor-product basis function coefficients
	coefficients = 0

	# logarithmic bias
	bias = 0

	# extent of supported region
	extents = []

	# geometry type
	geometry = 2
