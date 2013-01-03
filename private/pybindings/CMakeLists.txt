
if(NUMPY_FOUND)
i3_add_pybindings(photospline
	I3SplineTable.cxx
	module.cxx
	USE_TOOLS boost python numpy 
	USE_PROJECTS icetray photospline
)
endif(NUMPY_FOUND)

