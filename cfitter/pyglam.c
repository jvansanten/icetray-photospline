#include <Python.h>
#include <Numeric/arrayobject.h>
#include <stdio.h>

#include "glam.h"

/* Sparse GLAM Python interface */

static PyObject *pybox(PyObject *self, PyObject *args);

static PyMethodDef methods[] = {
	{ "box", pybox, METH_VARARGS },
	{ NULL, NULL }
};

void initspglam()
{
	import_array();
	Py_InitModule("spglam", methods);
}

static cholmod_sparse *
numpy2d_to_sparse(PyArrayObject *a, cholmod_common *c)
{
	cholmod_dense ad;

	ad.nrow = a->dimensions[0];
	ad.ncol = a->dimensions[1];
	ad.nzmax = ad.nrow * ad.ncol;
	ad.d = ad.nrow;
	ad.x = a->data;
	ad.z = NULL;
	ad.xtype = CHOLMOD_REAL;
	ad.dtype = CHOLMOD_DOUBLE;

	return cholmod_dense_to_sparse(&ad, 1, c);
}

static PyArrayObject *
numpy_sparse_to_2d(cholmod_sparse *a, cholmod_common *c)
{
	int dimensions[2];
	PyArrayObject *out;
	cholmod_dense *ad;

	dimensions[0] = a->nrow;
	dimensions[1] = a->ncol;
	out = (PyArrayObject *)PyArray_FromDims(2, dimensions,
	    PyArray_DOUBLE);

	ad = cholmod_sparse_to_dense(a, c);
	memcpy(out->data, ad->x, sizeof(double) * ad->nrow * ad->ncol);
	cholmod_free_dense(&ad, c);

	return out;
}
	

static PyObject *pybox(PyObject *self, PyObject *args)
{
	PyObject *xa, *xb;
	PyArrayObject *a, *b, *result_array;
	cholmod_common c;
	cholmod_sparse *am, *bm, *result;

	if (!PyArg_ParseTuple(args, "OO", &xa, &xb))
		return NULL;

	a = (PyArrayObject *)PyArray_ContiguousFromObject(xa,
	    PyArray_DOUBLE, 2, 2);
	b = (PyArrayObject *)PyArray_ContiguousFromObject(xb,
	    PyArray_DOUBLE, 2, 2);

	if (a == NULL || b == NULL) {
		PyErr_SetString(PyExc_ValueError,
		    "arrays must be two-dimensional and of type float");
		return NULL;
	}

	cholmod_start(&c);

	am = numpy2d_to_sparse(a, &c);
	bm = numpy2d_to_sparse(b, &c);

	Py_DECREF(a);
	Py_DECREF(b);

	result = box(am, bm, &c);
	cholmod_free_sparse(&am, &c);
	cholmod_free_sparse(&bm, &c);

	result_array = numpy_sparse_to_2d(result, &c);

	cholmod_free_sparse(&result, &c);
	cholmod_finish(&c);

	return PyArray_Return(result_array);
}
