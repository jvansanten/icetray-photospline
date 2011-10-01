
#define NO_IMPORT_ARRAY /* Just use the headers */
#define PY_ARRAY_UNIQUE_SYMBOL photospline_PyArray_API
#include <numpy/ndarrayobject.h>

#include "photospline/I3SplineTable.h"

namespace bp = boost::python;

#define PY_TYPESTRING(pyobj) \
	pyobj.ptr()->ob_type->tp_name

static double
splinetableeval(I3SplineTable &self, bp::object coordinates)
{
	double retvalue(NAN);
	double *coord_ptr;

	PyObject *coords;	
	coords = PyArray_ContiguousFromObject(coordinates.ptr(), NPY_DOUBLE, 0, 0);
	if (!coords) {
		PyErr_Format(PyExc_ValueError, "Can't convert object of type"
		    "'%s' to an array of doubles!", PY_TYPESTRING(coordinates));
		bp::throw_error_already_set();
	}
	coord_ptr = (double*)PyArray_DATA(coords);
	
	self.Eval(coord_ptr, &retvalue);

	return retvalue;
}

void register_I3SplineTable() {
	bp::class_<I3SplineTable, boost::shared_ptr<I3SplineTable> >
	    ("I3SplineTable", bp::init<const std::string&>(bp::arg("path")))
	    .def("eval", splinetableeval, bp::args("coordinates"), "Evaluate "
	        "the spline surface at the given coordinates.")
	;
}

