
if(NUMPY_FOUND)
i3_add_pybindings(photospline
	#cxx/I3SplineTable.cxx
	I3SplineTable.cxx
	module.cxx
	USE_TOOLS boost python numpy log4cplus
	USE_PROJECTS icetray photospline
)
endif(NUMPY_FOUND)

