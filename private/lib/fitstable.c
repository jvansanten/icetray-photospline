#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <fitsio.h>
#include <fitsio2.h>

#include "photospline/splinetable.h"

static int parsefitstable(fitsfile *fits, struct splinetable *table);

int
readsplinefitstable(const char *path, struct splinetable *table)
{
	fitsfile *fits;
	int error = 0;

	memset(table, 0, sizeof(struct splinetable));

	fits_open_file(&fits, path, READONLY, &error);
	if (error != 0)
		return (error);

	error = parsefitstable(fits, table);
	fits_close_file(fits, &error);
	fits_report_error(stderr, error);

	return (error);
}

static void
splinetable_free_aux(struct splinetable *table)
{
	int i;

	for (i = 0; i < table->naux; i++) {
		free(table->aux[i][0]);
		free(table->aux[i][1]);
		free(table->aux[i]);
	}

	free(table->aux);
	table->naux = 0;
	table->aux = NULL;
}

void splinetable_free(struct splinetable *table)
{
	int i;

	splinetable_free_aux(table);
	if (table->knots) {
		for (i = 0; i < table->ndim; i++)
			free(table->knots[i]);
		free(table->knots);
	}
	free(table->nknots);
	free(table->naxes);
	free(table->coefficients);
	free(table->periods);
	free(table->strides);
	
	if (table->extents) {
		free(table->extents[0]);
		free(table->extents);
	}
}

static int
parsefitstable(fitsfile *fits, struct splinetable *table)
{
	int error = 0;
	int hdus, type, i, nkeys;
	size_t arraysize;
	long *fpixel;

	fits_get_num_hdus(fits, &hdus, &error);
	fits_movabs_hdu(fits, 1, &type, &error);
	if (error != 0)
		return (error);

	if (type != IMAGE_HDU)
		return (ENOENT);

	/*
	 * Read header information
	 */

	fits_get_img_dim(fits, &table->ndim, &error);
	if (error != 0)
		return (error);

	/*
	 * Read in any auxiliary keywords.
	 */
	nkeys = 0;
	fits_get_hdrspace(fits, &nkeys, NULL, &error);
	if (nkeys > 0) {
		char key[FLEN_KEYWORD], value[FLEN_VALUE];
		int keylen, valuelen;
		table->aux = calloc(sizeof(char**), nkeys);
		i = 0;
		int j = 1;
		for ( ; (i < nkeys) && (j-1 < nkeys); j++) {
			error = 0;
			fits_read_keyn(fits, j, key, value, NULL, &error);
			if (error != 0)
				continue;
			if (strncmp("TYPE", key, 4) == 0 ||
			    strncmp("ORDER", key, 5) == 0 || 
			    strncmp("NAXIS", key, 5) == 0 ||
			    strncmp("PERIOD", key, 6) == 0 ||
			    strncmp("EXTEND", key, 6) == 0)
				continue;

			keylen = strlen(key) + 1;
			valuelen = strlen(value) + 1;
			table->aux[i] = calloc(sizeof(char*), 2);
			table->aux[i][0] = calloc(sizeof(char), keylen);
			table->aux[i][1] = calloc(sizeof(char), valuelen);
			memcpy(table->aux[i][0], key, keylen);
			memcpy(table->aux[i][1], value, valuelen);
			i++;
		}
		table->aux = realloc(table->aux, i*sizeof(char**));
		table->naux = i;
	} else {
		table->aux = NULL;
		table->naux = 0;
	}

	table->order = malloc(sizeof(table->order[i])*table->ndim);
	fits_read_key(fits, TINT, "ORDER", &table->order[0], NULL, &error);
	if (error != 0) {
		error = 0;
		
		for (i = 0; i < table->ndim; i++) {
			char name[255];
			sprintf(name,"ORDER%d",i);
			fits_read_key(fits, TINT, name, &table->order[i],
			    NULL, &error);
		}
	} else {
		for (i = 1; i < table->ndim; i++)
			table->order[i] = table->order[0];
	}

	if (error != 0)
		return (error);

	error = 0;
	table->periods = malloc(sizeof(table->periods[i])*table->ndim);
	for (i = 0; i < table->ndim; i++) {
		char name[255];
		sprintf(name,"PERIOD%d",i);
		fits_read_key(fits, TDOUBLE, name, &table->periods[i], NULL, &error);
		/*
		 * If the PERIOD keys cannot be read, just interpret this as
		 * a non-periodic table.
		 */
		if (error != 0) {
			table->periods[i] = 0;
			error = 0;
		}
	}

	/*
	 * Read the coefficient table
	 */

	table->naxes = malloc(sizeof(long)*table->ndim);
	fits_get_img_size(fits, table->ndim, table->naxes, &error);

	/*
	 * FITS multidimensional arrays are stored as FORTRAN arrays,
	 * not C arrays, so we need to swizzle the matrix into being
	 * a C array. Or we should. Instead, PyFITS, which writes these
	 * files, writes a C array, but with the axis counts transposed.
	 * Fix it.
	 */
	{
		long *naxestmp = malloc(sizeof(long)*table->ndim);
		for (i = 0; i < table->ndim; i++)
			naxestmp[i] = table->naxes[table->ndim - i - 1];

		free(table->naxes);
		table->naxes = naxestmp;
	}

	/* Compute the total array size and the strides into each dimension */
	table->strides = malloc(sizeof(unsigned long)*table->ndim);
	table->strides[table->ndim - 1] = arraysize = 1;
	for (i = table->ndim-1; i >= 0; i--) {
		arraysize *= table->naxes[i];
		if (i > 0)
			table->strides[i-1] = arraysize;
	}
	table->coefficients = malloc(sizeof(float)*arraysize);

	fpixel = malloc(sizeof(long)*table->ndim);
	for (i = 0; i < table->ndim; i++)
		fpixel[i] = 1;

	fits_read_pix(fits, TFLOAT, fpixel, arraysize, NULL,
	    table->coefficients, NULL, &error);

	free(fpixel);

	if (error != 0) {
		fprintf(stderr, "Error reading table coefficients\n");
		splinetable_free(table);
		return (error);
	}

	/*
	 * Read the knot vectors, which are stored one each in extension
	 * HDUs
	 */

	table->knots = malloc(sizeof(table->knots[0])*table->ndim);
	table->nknots = malloc(sizeof(table->nknots[0])*table->ndim);

	for (i = 0; i < table->ndim; i++) {
		char hduname[255];
		long fpix = 1;
		sprintf(hduname,"KNOTS%d",i);

		fits_movnam_hdu(fits, IMAGE_HDU, hduname, 0, &error);
		fits_get_img_size(fits, 1, &table->nknots[i], &error);
		if (error != 0) {
			fprintf(stderr, "Error reading knot vector %d\n", i);
			break;
		}

		table->knots[i] = malloc(sizeof(double)*table->nknots[i]);
		fits_read_pix(fits, TDOUBLE, &fpix, table->nknots[i], NULL,
		    table->knots[i], NULL, &error);
	}

	/*
	 * Read the axes extents, stored in a single extension HDU.
	 */

	table->extents = malloc(sizeof(double*) * table->ndim);
	table->extents[0] = malloc(sizeof(double) * 2 * table->ndim);

	for (i = 1; i < table->ndim; i++) {
		table->extents[i] = &table->extents[0][2*i];
	}

	long n_extents = 0;
	long fpix = 1;
	int ext_error = 0;
	fits_movnam_hdu(fits, IMAGE_HDU, "EXTENTS", 0, &ext_error);
	fits_get_img_size(fits, 1, &n_extents, &ext_error);
	if (n_extents != 2*table->ndim)
		ext_error = 1;

	if (ext_error != 0) { /* No extents. Make up some reasonable ones. */
		for (i = 0; i < table->ndim; i++) {
			table->extents[i][0] = table->knots[i][table->order[i]];
			table->extents[i][1] = table->knots[i][table->nknots[i]
			    - table->order[i] - 1];
		}
	} else {
		fits_read_pix(fits, TDOUBLE, &fpix, n_extents, NULL,
		     table->extents[0], NULL, &ext_error);
	}

	return (error);
}

char *
splinetable_get_key(struct splinetable *table, const char *key)
{
	int i = 0;
	char *value = NULL;
	
	for ( ; i < table->naux; i++) {
		if (strcmp(key, table->aux[i][0]) == 0) {
			value = table->aux[i][1];
		}
	}

	return (value);
}

int
splinetable_read_key(struct splinetable *table, splinetable_dtype type,
    const char *key, void *result)
{
	int error = 0;
	char *value = splinetable_get_key(table, key);

	if (!value)
		return (-1);

	switch (type) {
		case SPLINETABLE_INT:
			ffc2i(value, (long*)result, &error);
			break;
		case SPLINETABLE_DOUBLE:
			ffc2d(value, (double*)result, &error);
			break;
		default:
			error = BAD_DATATYPE;
	}

	if (error != 0)
		return (-1);
	else
		return (0);
	
}
