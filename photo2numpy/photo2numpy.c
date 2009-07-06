#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <stdio.h>
#include <limits.h>
#include <photonics.h>
#include <level2_reader.h>

#define L1_MAXDIM 6
#define L2_MAXDIM 4
#define TABLE_CHUNKSIZE 1024

static PyObject *photol1_chunks_to_numpy(Header_type *photoheader, FILE *table);
static PyObject *photol2_chunks_to_numpy(Level2_header_type *photoheader,
    FILE *table);

#define min(x,y) (((x) < (y)) ? (x) : (y))

/* Raw Photonics table Python interface */

static PyObject *readl1table(PyObject *self, PyObject *args);
static PyObject *readl2table(PyObject *self, PyObject *args);

static PyMethodDef methods[] = {
	{ "readl1", readl1table, METH_VARARGS },
	{ "readl2", readl2table, METH_VARARGS },
	{ NULL, NULL }
};

void initphoto2numpy()
{
	import_array();
	Py_InitModule("photo2numpy", methods);
}

static PyObject *readl1table(PyObject *self, PyObject *args)
{
	const char *path;
	PyObject *main_array, *stats_array;
	PyObject *coords[L1_MAXDIM], *coords_tuple;
	PyObject *binwidths[L1_MAXDIM], *binwidths_tuple;
	PyObject *result;
	char coordstr[L1_MAXDIM];
	FILE *table;
	int i, j, ndim;

	/* Glue variables for Photonics */
	Io_type io={0,{0,0,0,0,0,0},0,NULL};
	Header_type photoheader;
	Geo_type geo;

	result = NULL;

	if (!PyArg_ParseTuple(args, "s", &path))
		return NULL;

	printf("Opening photon table: %s\n",path);

	/* This mess is copied from the poorly named photonics function
	 * test_input() */

	table = fopen(path, "rb");
	if (table == NULL) {
		PyErr_SetString(PyExc_IOError,
		    "opening table failed");
		goto exit;
	}

	io.h = malloc(sizeof(Header_type));
	if (!read_header(io.h, table)) {
		PyErr_SetString(PyExc_ValueError,
		    "parsing table header failed");
		goto exit;
	}

	#ifdef SUPPORT_BIGENDIAN
		if(!isLittleEndian() && checkMetaHeadLittle(&(io.h->MetaHead)))
			byteswap32_l1header(io.h);
	#endif

	io.h_offset = sizeof(Header_type);
	io.offset[VARS-1] = sizeof(float);
	for(i = VARS-1; i > 0; --i)
		io.offset[i-1] = io.offset[i]*io.h->n[i];
	io.n_chunk = io.h->n[VARS-1];

	/* The following magic is from photo2ascii */

	copy_header(io.h,&photoheader);

	for (i = 0; i < VARS; i++) {
		photoheader.maxes[i][0] = photoheader.limits[i][0];
		photoheader.maxes[i][1] = photoheader.limits[i][1];
	}

	if (!set_up_geometry(0,&photoheader,&geo)) {
		PyErr_SetString(PyExc_ValueError,
		    "unable to parse photonics table geometry");
		goto exit;
	}

	/* Read out the main part of the table */

	main_array = (PyObject *)photol1_chunks_to_numpy(&photoheader, table);

	/* Check for statistics */
	if (photoheader.record_errors)
		stats_array = (PyObject *)photol1_chunks_to_numpy(
		    &photoheader, table);
	else
		stats_array = Py_None;

	/* Build up a tuple of the bin coordinates on each axis */
	ndim = 0;
	memset(coordstr, 0, sizeof(coordstr));
	for (i = 0; i < L1_MAXDIM; i++) {
		PyArrayObject *axis, *widths;
		npy_intp extent[1];

		if (photoheader.n[i] <= 1)
			continue;
		
		extent[0] = photoheader.n[i];
		axis = (PyArrayObject *)PyArray_SimpleNew(1, extent,
		    PyArray_DOUBLE);
		widths = (PyArrayObject *)PyArray_SimpleNew(1, extent,
		    PyArray_DOUBLE);

		for (j = 0; j < photoheader.n[i]; j++) {
			double val, width;

			switch (i) {
				case 5:
					val = (geo.timing[j] +
					    geo.timing[j+1])/2.;
					width = geo.timing[j+1] - geo.timing[j];
					break;
				case 1:
					if (photoheader.geo != CUBIC) {
						val = 180.0/M_PI *
						    (acos(geo.bracket[i][j]) +
						     acos(geo.bracket[i][j+1]));
						val /= 2.0;
						width = 180.0/M_PI *
						    (acos(geo.bracket[i][j+1]) -
						     acos(geo.bracket[i][j]));
						break;
					}

					/* In the cubic case, azimuth is just
					 * like other coords, so fallthrough. */
				default:
					val = (geo.bracket[i][j] +
					    geo.bracket[i][j+1])/2.;
					width = geo.bracket[i][j+1] -
					    geo.bracket[i][j];
			}

			((double *)(axis->data))[j] = val;
			((double *)(widths->data))[j] = width;
		}

		coords[ndim] = (PyObject *)axis;
		binwidths[ndim] = (PyObject *)widths;
		coordstr[ndim] = 'O';
		ndim++;
	}
	
	coords_tuple = Py_BuildValue(coordstr, coords[0], coords[1], coords[2],
	    coords[3], coords[4], coords[5]);
	binwidths_tuple = Py_BuildValue(coordstr, binwidths[0], binwidths[1],
	    binwidths[2], binwidths[3], binwidths[4], binwidths[5]);

	/* Now put together the final result */

	result = Py_BuildValue("OOOO", main_array, stats_array, coords_tuple,
	    binwidths_tuple);

    exit:
	free(io.h);
	if (table != NULL)
		fclose(table);

	return (result);
}

static PyObject *photol1_chunks_to_numpy(Header_type *photoheader,
    FILE *table)
{
	PyArrayObject *result;
	npy_intp dimensions[L1_MAXDIM];
	size_t valsleft;
	int ndim, i, j, valsread;
	float *array;

	/* Set up numpy array */
	ndim = 0;
	for (i = 0; i < L1_MAXDIM; i++) {
		/* Don't include dimensions that have only one element */
		if(photoheader->n[i] > 1) {
			dimensions[ndim] = photoheader->n[i];
			ndim++;
		}
	}

	result = (PyArrayObject *)PyArray_SimpleNew(ndim, dimensions,
	    PyArray_DOUBLE);

	/*
	 * Now loop through the file, copying TABLE_CHUNKSIZE values at a time.
	 */

	valsleft = 1;
	for (i = 0; i < L1_MAXDIM; i++)
		valsleft *= photoheader->n[i];

	i = 0;
	array = malloc(sizeof(float)*TABLE_CHUNKSIZE);
	errno = 0;
	while (valsleft > 0) {
		valsread = fread(array, sizeof(float),
		    min(TABLE_CHUNKSIZE, valsleft), table);
		if (valsread == 0 || ferror(table)) {
			if (valsread == 0)
				fprintf(stderr,"Unexpected end-of-file!\n");
			else
				fprintf(stderr,"Error while reading table: %s\n",
				    strerror(errno));
			Py_DECREF(result);
			result = NULL;
			goto exit;
		}

		for (j = 0; j < valsread; j++)
			((double *)(result->data))[i++] = array[j];

		valsleft -= valsread;
	}

	/* All done */
    exit:
	free(array);

	if (result == NULL)
		return (Py_None);

	return ((PyObject *)result);
}

static PyObject *readl2table(PyObject *self, PyObject *args)
{
	const char *path;
	PyObject *main_array, *stats_array;
	PyObject *coords[L2_MAXDIM], *coords_tuple;
	PyObject *binwidths[L2_MAXDIM], *binwidths_tuple;
	PyObject *result;
	char coordstr[L2_MAXDIM];
	FILE *table;
	int i, j, ndim;

	/* Glue variables for Photonics */
	Level2_header_type photoheader;
	Level2_geo_type geo;

	result = NULL;

	if (!PyArg_ParseTuple(args, "s", &path))
		return NULL;

	printf("Opening photon table: %s\n",path);

	/* This mess is copied from the poorly named photonics function
	 * test_input() */

	table = fopen(path, "rb");
	if (table == NULL) {
		PyErr_SetString(PyExc_IOError,
		    "opening table failed");
		goto exit;
	}

	if (fread(&photoheader, sizeof(Level2_header_type), 1, table) <= 0) {
		PyErr_SetString(PyExc_ValueError,
		    "reading table header failed");
		goto exit;
	}

	if (!level2_geometry(&photoheader, &geo)) {
		PyErr_SetString(PyExc_ValueError,
		    "parsing table header failed");
		goto exit;
	}

	/* Fudge the time axis for Level 2 ABS tables */
	if (photoheader.type == ABS)
		photoheader.n[3] = 1; /* Only 1 time bin */	

	/* Read out the main part of the table */

	main_array = (PyObject *)photol2_chunks_to_numpy(&photoheader, table);
	stats_array = Py_None;

	/* Build up a tuple of the bin coordinates on each axis */
	ndim = 0;
	memset(coordstr, 0, sizeof(coordstr));
	for (i = 0; i < L1_MAXDIM; i++) {
		PyArrayObject *axis, *widths;
		npy_intp extent[1];

		if (photoheader.n[i] <= 1)
			continue;
		
		extent[0] = photoheader.n[i];
		axis = (PyArrayObject *)PyArray_SimpleNew(1, extent,
		    PyArray_DOUBLE);
		widths = (PyArrayObject *)PyArray_SimpleNew(1, extent,
		    PyArray_DOUBLE);

		for (j = 0; j < photoheader.n[i]; j++) {
			double val, width;

			switch (i) {
				case 0:
					val = photoheader.range[0][0] +
					    (j + 0.5)*geo.d[0];
					width = geo.d[0];
					break;
				case 1:
					val = (j + 0.5)*geo.d[1];
					if (photoheader.d_scale == QUADRATIC) {
						val *= val;
						width = sqr((j + 1.5)*geo.d[1])
						    - val;
					} else {
						width = geo.d[1];
					}
					val += photoheader.range[1][0];
					break;
				case 2:
					val = photoheader.range[2][0] +
					    j * geo.d[2];
					width = geo.d[2];
					break;
				case 3:
					val = (j + 0.5)*geo.d[3];
					if (photoheader.t_scale == QUADRATIC) {
						val *= val;
						width = sqr((j + 1.5)*geo.d[3])
						    - val;
					} else {
						width = geo.d[3];
					}
					val += photoheader.range[1][0];
					break;
				default:
					PyErr_SetString(PyExc_ValueError,
					    "number of dimensions exceeds 4");
					goto exit;
			}

			((double *)(axis->data))[j] = val;
			((double *)(widths->data))[j] = width;
		}

		coords[ndim] = (PyObject *)axis;
		binwidths[ndim] = (PyObject *)widths;
		coordstr[ndim] = 'O';
		ndim++;
	}
	
	coords_tuple = Py_BuildValue(coordstr, coords[0], coords[1], coords[2],
	    coords[3], coords[4], coords[5]);
	binwidths_tuple = Py_BuildValue(coordstr, binwidths[0], binwidths[1],
	    binwidths[2], binwidths[3], binwidths[4], binwidths[5]);

	/* Now put together the final result */

	result = Py_BuildValue("OOOO", main_array, stats_array, coords_tuple,
	    binwidths_tuple);

    exit:
	if (table != NULL)
		fclose(table);

	return (result);
}

static PyObject *photol2_chunks_to_numpy(Level2_header_type *photoheader,
    FILE *table)
{
	PyArrayObject *result;
	npy_intp dimensions[L2_MAXDIM];
	size_t valsleft;
	int ndim, i, j, valsread;
	float *array;

	/* Set up numpy array */
	ndim = 0;
	for (i = 0; i < L2_MAXDIM; i++) {
		/* Don't include dimensions that have only one element */
		if(photoheader->n[i] > 1) {
			dimensions[ndim] = photoheader->n[i];
			ndim++;
		}
	}

	result = (PyArrayObject *)PyArray_SimpleNew(ndim, dimensions,
	    PyArray_DOUBLE);

	/*
	 * Now loop through the file, copying TABLE_CHUNKSIZE values at a time.
	 */

	valsleft = 1;
	for (i = 0; i < L2_MAXDIM; i++)
		valsleft *= photoheader->n[i];

	i = 0;
	array = malloc(sizeof(float)*TABLE_CHUNKSIZE);
	errno = 0;
	while (valsleft > 0) {
		valsread = fread(array, sizeof(float),
		    min(TABLE_CHUNKSIZE, valsleft), table);
		if (valsread == 0 || ferror(table)) {
			if (valsread == 0)
				fprintf(stderr,"Unexpected end-of-file!\n");
			else
				fprintf(stderr,"Error while reading table: %s\n",
				    strerror(errno));
			Py_DECREF(result);
			result = NULL;
			goto exit;
		}

		for (j = 0; j < valsread; j++)
			((double *)(result->data))[i++] = array[j];

		valsleft -= valsread;
	}

	/* All done */
    exit:
	free(array);

	if (result == NULL)
		return (Py_None);

	return ((PyObject *)result);
}



