
#define PY_ARRAY_UNIQUE_SYMBOL photospline_PyArray_API

#include <icetray/load_project.h>
#include <numpy/ndarrayobject.h>

void register_I3SplineTable();

#if PY_VERSION_HEX >= 0x03000000
static PyObject *hack_numpy()
#else
static void hack_numpy()
#endif
{
	import_array();
#if PY_VERSION_HEX >= 0x03000000
	return NULL;
#endif
}

I3_PYTHON_MODULE(photospline)
{
	load_project("photospline", false);
	register_I3SplineTable();

	// boost::python::import("numpy");
	hack_numpy();
}
