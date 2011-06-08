#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <stdio.h>
#include <limits.h>
#include <photonics.h>
#ifndef SKIP_LEVEL2
#include <level2_reader.h>
#endif
#include "hobo_mmap.h"

#define L1_MAXDIM 6
#define L2_MAXDIM 4
#define TABLE_CHUNKSIZE 1024

static PyObject *photol1_chunks_to_numpy(Header_type *photoheader, FILE *table);
static PyObject *photol1_map_to_numpy(Header_type *photoheader, FILE *table);
#ifndef SKIP_LEVEL2
static PyObject *photol2_chunks_to_numpy(Level2_header_type *photoheader,
    FILE *table);
#endif

#define min(x,y) (((x) < (y)) ? (x) : (y))

/* Raw Photonics table Python interface */

static PyObject *readl1table(PyObject *self, PyObject *args);
#ifndef SKIP_LEVEL2
static PyObject *readl2table(PyObject *self, PyObject *args);
#endif

static PyMethodDef methods[] = {
	{ "readl1", readl1table, METH_VARARGS },
#ifndef SKIP_LEVEL2
	{ "readl2", readl2table, METH_VARARGS },
#endif
	{ NULL, NULL }
};

void initphoto2numpy(void)
{
	PyType_Ready(&mmap_wrapper_type);

	import_array();
	Py_InitModule("photo2numpy", methods);
}

static PyObject *readl1table(PyObject *self, PyObject *args)
{
	const char *path;
	PyObject *main_array, *stats_array;
	PyObject *coords[L1_MAXDIM], *coords_tuple;
	PyObject *binwidths[L1_MAXDIM], *binwidths_tuple;
	PyObject *header_dict;
	PyObject *result;
	PyObject *py_do_map;
	char coordstr[L1_MAXDIM+1];
	FILE *table;
	int i, j, ndim, do_map;

	/* Glue variables for Photonics */
	Io_type io={0,{0,0,0,0,0,0},0,NULL};
	Header_type photoheader;
	Geo_type geo;

	result = NULL;
	py_do_map = Py_False;

	if (!PyArg_ParseTuple(args, "sO:readl1", &path, &py_do_map))
		return NULL;

	do_map = PyObject_IsTrue(py_do_map);

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

	header_dict = PyDict_New();
	PyDict_SetItemString(header_dict, "n_photon",
	    PyInt_FromLong(io.h->n_photon));
	PyDict_SetItemString(header_dict, "efficiency",
	    PyInt_FromLong((long)(io.h->efficiency)));
	PyDict_SetItemString(header_dict, "geometry",
	    PyInt_FromLong((long)(io.h->geo)));
	PyDict_SetItemString(header_dict, "zenith",
	    PyFloat_FromDouble((double)(io.h->angle)));
	PyDict_SetItemString(header_dict, "z",
	    PyFloat_FromDouble((double)(io.h->depth)));

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

	/* 
	 * After read_header(), the file handle is positioned at the head
	 * of the main data table. Read out the main part of the table.
	 */

	if (do_map)
		main_array = photol1_map_to_numpy(&photoheader, table);
	else
		main_array = photol1_chunks_to_numpy(&photoheader, table);

	/* 
	 * If the table was generated with statistics, there is a second data
	 * table following the main table and containing either the number of
	 * photons tracked through the bin or the sum of the squares of the
	 * weights assigned to those photons. Reading the main table sets the
	 * file position to the head of this second table, so we simply turn the
	 * crank again.
	 */
	if (photoheader.record_errors) {
		if (do_map)
			stats_array = photol1_map_to_numpy(&photoheader, table);
		else
			stats_array = photol1_chunks_to_numpy(&photoheader,
			    table);
	} else {
		stats_array = Py_None;
	}

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

	result = Py_BuildValue("OOOOO", main_array, stats_array, coords_tuple,
	    binwidths_tuple, header_dict);

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

static PyObject *photol1_map_to_numpy(Header_type *photoheader,
    FILE *table)
{
	PyArrayObject *result;
	PyObject *mmap;
	npy_intp dimensions[L1_MAXDIM];
	size_t table_size;
	off_t pos;
	int i, ndim;
	void *data;

	/* Set up numpy array */
	ndim = 0;
	table_size = 1;
	for (i = 0; i < L1_MAXDIM; i++) {
		/* Don't include dimensions that have only one element */
		if(photoheader->n[i] > 1) {
			dimensions[ndim] = photoheader->n[i];
			table_size *= photoheader->n[i];
			ndim++;
		}
	}

	result = NULL;
	/* The file handle is positioned at the head of the table to be read */
	pos = ftello(table);
	
	/* We know where the table starts; seek to the end to emulate fread() */
	fseeko(table, table_size * sizeof(float), SEEK_CUR);

	mmap = (PyObject*)mmap_wrapper_new(fileno(table), pos,
	    table_size * sizeof(float), &data);

	if (!mmap) goto exit;

	result = (PyArrayObject*)PyArray_SimpleNewFromData(ndim, dimensions,
	    PyArray_FLOAT, data);

	if (!result) goto exit;

	Py_INCREF(mmap);
	result->base = mmap;

	/* All done */
    exit:

	if (result == NULL)
		return (Py_None);

	return ((PyObject *)result);
}

#ifndef SKIP_LEVEL2
static PyObject *readl2table(PyObject *self, PyObject *args)
{
	const char *path;
	PyObject *main_array, *stats_array;
	PyObject *coords[L2_MAXDIM], *coords_tuple;
	PyObject *binwidths[L2_MAXDIM], *binwidths_tuple;
	PyObject *header_dict;
	PyObject *result;
	char coordstr[L2_MAXDIM+1];
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
	
	header_dict = PyDict_New();
#if 0
	/* No useful normalization info */
	PyDict_SetItemString(header_dict, "n_photon",
	    PyInt_FromLong(io.h->n_photon));
#endif
	PyDict_SetItemString(header_dict, "efficiency",
	    PyInt_FromLong((long)(photoheader.efficiency)));
	/* NB: Level2 tables are always cylindrical */
	PyDict_SetItemString(header_dict, "geometry",
	    PyInt_FromLong(CYLINDRICAL));
	PyDict_SetItemString(header_dict, "zenith",
	    PyFloat_FromDouble((double)(photoheader.theta)));
	PyDict_SetItemString(header_dict, "z",
	    PyFloat_FromDouble((double)(photoheader.z0)));

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

	result = Py_BuildValue("OOOOO", main_array, stats_array, coords_tuple,
	    binwidths_tuple, header_dict);

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
#endif



