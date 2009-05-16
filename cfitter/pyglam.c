#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <stdio.h>
#include <limits.h>

#include "glam.h"

/* Sparse GLAM Python interface */

static PyObject *pybox(PyObject *self, PyObject *args);
static PyObject *pyrho(PyObject *self, PyObject *args);
static PyObject *pyfit(PyObject *self, PyObject *args);

static PyMethodDef methods[] = {
	{ "box", pybox, METH_VARARGS },
	{ "rho", pyrho, METH_VARARGS },
	{ "fit", pyfit, METH_VARARGS },
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

	sp = cholmod_l_dense_to_sparse(&ad, 1, c);
	
	/* Correct for row-major/column-major ordering issues */
	spt = cholmod_l_transpose(sp, 1, c);
	cholmod_l_free_sparse(&sp, c);

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

	at = cholmod_l_transpose(a, 1, c); /* Fix row-major/column-major */
	ad = cholmod_l_sparse_to_dense(at, c);
	cholmod_l_free_sparse(&at, c);
	memcpy(out->data, ad->x, sizeof(double) * ad->nrow * ad->ncol);
	cholmod_l_free_dense(&ad, c);

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

	cholmod_l_start(&c);

	am = numpy2d_to_sparse(a, &c);
	bm = numpy2d_to_sparse(b, &c);

	Py_DECREF(a);
	Py_DECREF(b);

	result = box(am, bm, &c);
	cholmod_l_free_sparse(&am, &c);
	cholmod_l_free_sparse(&bm, &c);

	if (result == NULL) {
		PyErr_SetString(PyExc_ValueError,
		    "arrays must have compatible dimensions");
		cholmod_l_finish(&c);
		return NULL;
	}

	result_array = numpy_sparse_to_2d(result, &c);

	cholmod_l_free_sparse(&result, &c);
	cholmod_l_finish(&c);

	return PyArray_Return(result_array);
}

void print_ndsparse_py(struct ndsparse *a) {
	PyObject *py = (PyObject *)numpy_ndsparse_to_ndarray(a);
	PyObject_Print(py, stdout, 0);
	Py_DECREF(py);
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

	cholmod_l_start(&c);

	am = numpy2d_to_sparse(a, &c);
	err = slicemultiply(&nd, am, dim, &c);

	cholmod_l_free_sparse(&am, &c);
	cholmod_l_finish(&c);

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

static PyObject *pyfit(PyObject *self, PyObject *args)
{
	PyObject *z, *w, *coords, *knots, *periods;
	PyArrayObject *result;
	PyArrayObject *z_arr;
	struct ndsparse data;
	struct splinetable out;
	double *data_arr, *weights;
	double **c_coords;
	cholmod_common c;
	int order;
	double smooth;
	int i, j, k, elements, err;
	int *moduli;

	/* Initialize a few things to NULL */
	moduli = NULL;
	weights = NULL;
	z_arr = NULL;
	result = NULL;
	memset(&out, 0, sizeof(out));

	/* Parse our arguments from Python land */
	if (!PyArg_ParseTuple(args, "OOOOidO", &z, &w, &coords, &knots,
	    &order, &smooth, &periods))
		return NULL;

	/* Parse weights first to avoid storing data with 0 weight */
	err = numpynd_to_ndsparse(w, &data);
	if (err == 0) {
		z_arr = (PyArrayObject *)PyArray_ContiguousFromObject(z,
		    PyArray_DOUBLE, data.ndim, data.ndim);
	}

	if (err != 0 || z_arr == NULL) {
		PyErr_SetString(PyExc_ValueError,
		    "could not decode arrays");
		return NULL;
	}

	/* Now validate the input */

	result = NULL;
	for (i = 0; i < data.ndim; i++) {
		if (z_arr->dimensions[i] != data.ranges[i]) {
			Py_DECREF(z_arr);
			PyErr_SetString(PyExc_ValueError,
			    "weight and data array dimensions do not match");
			goto exit;
		}
	}

	/* Set up the data array */
	moduli = malloc(sizeof(int) * data.ndim);
	moduli[data.ndim-1] = 1;
	for (i = data.ndim-2; i >= 0; i--)
		moduli[i] = moduli[i+1]*data.ranges[i+1];

	data_arr = calloc(data.rows,sizeof(double));
	for (i = 0; i < data.rows; i++) {
		k = 0;
		for (j = 0; j < data.ndim; j++)
			k += moduli[j]*data.i[j][i];
		data_arr[i] = ((double *)(z_arr->data))[k];
	}
	free(moduli);

	/* Swap the data into the ndsparse structure to satisfy glamfit's
	 * calling convention */
	weights = data.x;
	data.x = data_arr;

	/* We don't need the Python structure anymore */
	Py_DECREF(z_arr);

	/* Check knot and coords for consistency */
	if (!PySequence_Check(knots) || data.ndim != PySequence_Length(knots)) {
		PyErr_SetString(PyExc_TypeError,
		    "knots must be a sequence with one row for each dimension");
		goto exit;
	}
	if (!PySequence_Check(knots) || data.ndim != PySequence_Length(knots)) {
		PyErr_SetString(PyExc_TypeError,
		    "coord must be a sequence with one row for each dimension");
		goto exit;
	}

	/* Start setting up the spline table */
	out.ndim = data.ndim;
	out.order = order;

	out.knots = calloc(out.ndim,sizeof(double *));
	out.nknots = calloc(out.ndim,sizeof(long));
	for (i = 0; i < PySequence_Length(knots); i++) {
		PyArrayObject *knot_vec;
		knot_vec = (PyArrayObject *)PyArray_ContiguousFromObject(
		    PySequence_GetItem(knots, i),
		    PyArray_DOUBLE, 1, 1);

		if (knot_vec == NULL) {
			PyErr_SetString(PyExc_TypeError,
			    "knots cannot be read as arrays");
			goto exit;
		}

		out.nknots[i] = knot_vec->dimensions[0];
		out.knots[i] = calloc(out.nknots[i], sizeof(double));
		memcpy(out.knots[i], knot_vec->data,
		    out.nknots[i] * sizeof(double));

		Py_DECREF(knot_vec);
	}
	
	c_coords = malloc(sizeof(double *)*out.ndim);
	for (i = 0; i < PySequence_Length(coords); i++) {
		PyArrayObject *coord_vec;
		coord_vec = (PyArrayObject *)PyArray_ContiguousFromObject(
		    PySequence_GetItem(coords, i),
		    PyArray_DOUBLE, 1, 1);

		if (coord_vec == NULL) {
			PyErr_SetString(PyExc_TypeError,
			    "coords cannot be read as arrays");
			goto exit;
		}

		if (coord_vec->dimensions[0] != data.ranges[i]) {
			PyErr_SetString(PyExc_ValueError,
			    "wrong number of coords");
			Py_DECREF(coord_vec);
			goto exit;
		}

		c_coords[i] = calloc(data.ranges[i], sizeof(double));
		memcpy(c_coords[i], coord_vec->data,
		    data.ranges[i] * sizeof(double));

		Py_DECREF(coord_vec);
	}

	/* Do the fit */
	cholmod_l_start(&c);
	glamfit(&data, weights, c_coords, &out, smooth, order, 1, &c);
	cholmod_l_finish(&c);

	/* Now process the splinetable into a numpy array */
	elements = 1;
	for (i = 0; i < out.ndim; i++)
		elements *= out.naxes[i];
	result = (PyArrayObject *)PyArray_SimpleNew(out.ndim, out.naxes,
	    PyArray_DOUBLE);
	memcpy(result->data, out.coefficients, elements*sizeof(double));

   exit:
	for (i = 0; i < data.ndim; i++)
		free(data.i[i]);
	free(data.x); free(data.ranges);
	if (weights) free(weights);
	if (out.knots) {
		free(out.nknots);
		free(out.periods);
		for (i = 0; i < out.ndim; i++)
			free(out.knots[i]);
		free(out.knots);
	}
	if (out.coefficients) {
		free(out.naxes);
		free(out.coefficients);
	}

	return (PyObject *)result;
}


