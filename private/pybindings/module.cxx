
#define PY_ARRAY_UNIQUE_SYMBOL photospline_PyArray_API

#include <icetray/load_project.h>
#include <numpy/ndarrayobject.h>

void register_I3SplineTable();

I3_PYTHON_MODULE(photospline)
{
	load_project("photospline", false);
	// boost::python::import("numpy");
	import_array();
	
	register_I3SplineTable();
}
