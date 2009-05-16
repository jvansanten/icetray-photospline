#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <stdio.h>
#include <limits.h>

#include "glam.h"

/* Sparse GLAM Python interface */

static PyObject *pybox(PyObject *self, PyObject *args);
static PyObject *pyrho(PyObject *self, PyObject *args);

static PyMethodDef methods[] = {
	{ "box", pybox, METH_VARARGS },
	{ "rho", pyrho, METH_VARARGS },
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
	cholmod_sparse *sp, *spt;

	ad.nrow = a->dimensions[1];
	ad.ncol = a->dimensions[0];
	ad.nzmax = ad.nrow * ad.ncol;
	ad.d = ad.nrow;
	ad.x = a->data;
	ad.z = NULL;
	ad.xtype = CHOLMOD_REAL;
	ad.dtype = CHOLMOD_DOUBLE;

	sp = cholmod_dense_to_sparse(&ad, 1, c);
	
	/* Correct for row-major/column-major ordering issues */
	spt = cholmod_transpose(sp, 1, c);
	cholmod_free_sparse(&sp, c);

	return spt;
}

static PyArrayObject *
numpy_sparse_to_2d(cholmod_sparse *a, cholmod_common *c)
{
	npy_intp dimensions[2];
	PyArrayObject *out;
	cholmod_dense *ad;
	cholmod_sparse *at;

	dimensions[0] = a->nrow;
	dimensions[1] = a->ncol;
	out = (PyArrayObject *)PyArray_SimpleNew(2, dimensions,
	    PyArray_DOUBLE);

	at = cholmod_transpose(a, 1, c); /* Fix row-major/column-major */
	ad = cholmod_sparse_to_dense(at, c);
	cholmod_free_sparse(&at, c);
	memcpy(out->data, ad->x, sizeof(double) * ad->nrow * ad->ncol);
	cholmod_free_dense(&ad, c);

	return out;
}

static int
numpynd_to_ndsparse(PyObject *in, struct ndsparse *out)
{
	PyArrayObject *inar;
	int i, j, elements, coord, currow;
	int *moduli;

	inar = (PyArrayObject *)PyArray_ContiguousFromObject(in,
	    PyArray_DOUBLE, 1, INT_MAX);
	if (inar == NULL)
		return -1;

	out->rows = 0;
	out->ndim = inar->nd;
	out->ranges = malloc(sizeof(int)*out->ndim);
	out->i = malloc(sizeof(int *)*out->ndim);
	moduli = malloc(sizeof(int)*out->ndim);

	elements = 1;
	for (i = 0; i < inar->nd; i++) {
		out->ranges[i] = inar->dimensions[i];
		elements *= inar->dimensions[i];
	}

	moduli[out->ndim-1] = 1;
	for (i = out->ndim-2; i >= 0; i--)
		moduli[i] = moduli[i+1]*out->ranges[i+1];
	
	/* Find how many non-zeros we have */
	for (i = 0; i < elements; i++) {
		if (((double *)(inar->data))[i] != 0.0)
			out->rows++;
	}

	/* Init coordinates */
	for (i = 0; i < inar->nd; i++)
		out->i[i] = malloc(sizeof(int)*out->rows);
	out->x = malloc(sizeof(double)*out->rows);

	currow = 0;
	for (i = 0; i < elements; i++)  {
		if (((double *)(inar->data))[i] == 0)
			continue;
		out->x[currow] = ((double *)(inar->data))[i];

		coord = i;
		for (j = 0; j < inar->nd; j++) {
			out->i[j][currow] = coord / moduli[j];
			coord = coord % moduli[j];
		}

		currow++;
	}
	Py_DECREF(inar);
	free(moduli);

	return 0;
}

static PyArrayObject *
numpy_ndsparse_to_ndarray(struct ndsparse *a)
{
	double *x;
	PyArrayObject *out;
	int moduli[a->ndim];
	npy_intp dimensions[a->ndim];
	int i, j, k, elements;

	/* Change the type of a->ranges to pacify numpy */
	for (i = 0; i < a->ndim; i++)
		dimensions[i] = a->ranges[i];

	out = (PyArrayObject *)PyArray_SimpleNew(a->ndim, dimensions,
	    PyArray_DOUBLE);

	moduli[a->ndim-1] = 1;
	for (i = a->ndim-2; i >= 0; i--)
		moduli[i] = moduli[i+1]*a->ranges[i+1];

	/* Find and initialize the data array to zero */
	x = (double *)out->data;
	elements = 1;
	for (i = 0; i < a->ndim; i++) 
		elements *= a->ranges[i];
	memset(x, 0, sizeof(double)*elements);

	for (i = 0; i < a->rows; i++) {
		k = 0;
		for (j = 0; j < a->ndim; j++)
			k += a->i[j][i]*moduli[j];
		
		x[k] = a->x[i];
	}

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

	if (result == NULL) {
		PyErr_SetString(PyExc_ValueError,
		    "arrays must have compatible dimensions");
		cholmod_finish(&c);
		return NULL;
	}

	result_array = numpy_sparse_to_2d(result, &c);

	cholmod_free_sparse(&result, &c);
	cholmod_finish(&c);

	return PyArray_Return(result_array);
}

void printndsparse(struct ndsparse *a) {
	double *x;
	int moduli[a->ndim];
	int i, j, k, elements;

	elements = 1;
	for (i = 0; i < a->ndim; i++) 
		elements *= a->ranges[i];
	moduli[a->ndim-1] = 1;
	for (i = a->ndim-2; i >= 0; i--)
		moduli[i] = moduli[i+1]*a->ranges[i+1];

	for (i = 0; i < a->ndim; i++)
		printf("Dimension %d: %d\n",i,a->ranges[i]);

	x = calloc(elements,sizeof(double));

	for (i = 0; i < a->rows; i++) {
		k = 0;
		for (j = 0; j < a->ndim; j++)
			k += a->i[j][i]*moduli[j];
		
		x[k] = a->x[i];
	}

	for (i = 0; i < elements; i++) {
		if (i % moduli[0] == 0 && i != 0) printf("\n");
		if (i % moduli[1] == 0 && i != 0) printf("\n");
		printf("%lf\t",x[i]);
	}
	printf("\n");

	free(x);
}

static PyObject *pyrho(PyObject *self, PyObject *args)
{
	/* This takes the reverse calling order in Python:
	 * (matrix, ndarray, dim) */

	PyObject *xa, *xb;
	PyArrayObject *a, *result_array;
	struct ndsparse nd;
	cholmod_sparse *am;
	cholmod_common c;
	int dim, i, err;

	if (!PyArg_ParseTuple(args, "OOi", &xa, &xb, &dim))
		return NULL;

	a = (PyArrayObject *)PyArray_ContiguousFromObject(xa,
	    PyArray_DOUBLE, 2, 2);
	if (a == NULL || numpynd_to_ndsparse(xb, &nd) != 0) {
		if (a != NULL) {
			Py_DECREF(a);
		}
	
		PyErr_SetString(PyExc_ValueError,
		    "could not decode arrays");
		return NULL;
	}

	cholmod_start(&c);

	am = numpy2d_to_sparse(a, &c);
	err = slicemultiply(&nd, am, dim, &c);

	cholmod_free_sparse(&am, &c);
	cholmod_finish(&c);

	if (err == 0) {
		result_array = numpy_ndsparse_to_ndarray(&nd);
	} else {
		PyErr_SetString(PyExc_ValueError,
		    "Dimensions do not match");
		result_array = NULL;
	}

	for (i = 0; i < nd.ndim; i++)
		free(nd.i[i]);
	free(nd.x); free(nd.ranges);

	return (PyObject *)result_array;
}

