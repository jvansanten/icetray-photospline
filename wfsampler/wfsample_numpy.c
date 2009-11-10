#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <stdio.h>
#include <limits.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "logsplinepdf.h"

#ifndef __FreeBSD__
void srandomdev() {
        unsigned long seed;
        int fd;

        fd = open("/dev/random",O_RDONLY);
        read(fd,&seed,sizeof(seed));
        srandom(seed);
        close(fd);
}
#endif

static PyObject *pandelsample_numpy(PyObject *self, PyObject *args);
static PyObject *wfsample_numpy(PyObject *self, PyObject *args, PyObject *kws);


static PyMethodDef methods[] = {
	{ "pandelsample", pandelsample_numpy, METH_VARARGS },
	{ "wfsample", (PyCFunction)wfsample_numpy, METH_KEYWORDS },
	{ NULL, NULL }
};

gsl_rng *rng = NULL;

void initwfsample()
{
	// Set up GSL RNG
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(gsl_rng_default);
	srandomdev();
	gsl_rng_set(rng,random());

	// Set up Python module
	import_array();
	Py_InitModule("wfsample", methods);
}

const double lambda = 71.; // meter
const double tau = 671.; // ns
const double x0 = 154.; // meter
const double c = 0.2998; // meter/ns
const double n = 1.34;

/*
 * GAMMA_PROPOSAL_SCALE can be used to widen the variance of the Pandel
 * about the mean, which makes sure the tails of the Pandel proposal
 * distribution are larger than the sampling distribution, preventing
 * underdispersion in the tails. Set to 1.0 by default, reasonable values
 * are between 1.0 and 1.5
 */
#define GAMMA_PROPOSAL_SCALE	1.0

/* THE FOLLOWING ARE NOT THREAD SAFE */
double distance;

static double pandel_sample()
{
	double gamma_scale = 1./(1./tau + c/n/x0);
	double gamma_shape = distance/lambda;
	
	/* Widen the distribution for better sampling */
	gamma_scale *= GAMMA_PROPOSAL_SCALE;
	gamma_shape /= GAMMA_PROPOSAL_SCALE;

	return gsl_ran_gamma(rng, gamma_shape, gamma_scale);
}

static double pandel_pdf(double t, double lastt)
{
	double gamma_scale = 1./(1./tau + c/n/x0);
	double gamma_shape = distance/lambda;

	/* Widen the distribution for better sampling */
	gamma_scale *= GAMMA_PROPOSAL_SCALE;
	gamma_shape /= GAMMA_PROPOSAL_SCALE;

	return gsl_ran_gamma_pdf(t, gamma_shape, gamma_scale);
}

static PyObject *pandelsample_numpy(PyObject *self, PyObject *args)
{
	int size, i;
	PyArrayObject *result;
	npy_intp extent[1];

	size = 200;
	if (!PyArg_ParseTuple(args, "di", &distance, &size))
		return NULL;

	extent[0] = size;
	result = (PyArrayObject *)PyArray_SimpleNew(1, extent, PyArray_DOUBLE);
	for (i = 0; i < size; i++)
		((double *)(result->data))[i] = pandel_sample();

	return (PyObject *)(result);
}

static PyObject *wfsample_numpy(PyObject *self, PyObject *args, PyObject *kws)
{
	const char *path;
	int size, burnin;
	PyArrayObject *result;
	PyObject *pycoords;
	struct splinetable table;
	npy_intp extent[1];
	double coords[6];

	char *kwlist[] = {"table", "coords", "size", "burnin", NULL};
	size = 200;
	burnin = 50;

	if (!PyArg_ParseTupleAndKeywords(args, kws, "sO|ii", kwlist, &path, 
	    &pycoords, &size, &burnin))
		return NULL;

	if (readsplinefitstable(path, &table) != 0) {
		PyErr_SetString(PyExc_IOError,
		    "opening table failed");
		return (NULL);
	}

	/* Parse the coordinates array */
	{
		PyArrayObject *py_arrcoords;
		int i,j;

		py_arrcoords = (PyArrayObject *)PyArray_ContiguousFromObject(
		    pycoords, PyArray_DOUBLE, 1, 1);
		if (py_arrcoords == NULL || py_arrcoords->dimensions[0] !=
		    table.ndim - 1) {
			PyErr_SetString(PyExc_ValueError,
			    "coordinate array must be one dimensional and have"
			    "one fewer value than the number of dimensions"
			    "in the table");
			return NULL;
		}

		coords[3] = 0;
		for (i = 0, j = 0; i < table.ndim; i++) {
			if (i == 3) continue;

			coords[i] = ((double *)py_arrcoords->data)[j++];
		}
		/* THIS ASSUMES CYLINDRICAL GEOMETRY */
		distance = hypot(coords[0], coords[2]);

		Py_DECREF((PyObject *)py_arrcoords);
	}

	extent[0] = size;
	result = (PyArrayObject *)PyArray_SimpleNew(1, extent, PyArray_DOUBLE);

	logsplinepdf_n_sample((double *)result->data, size, burnin, coords,
	    3 /* time */, &table, (1 << 3), &pandel_sample, &pandel_pdf, rng);

	return (PyObject *)(result);
}

