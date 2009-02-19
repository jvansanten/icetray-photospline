#include <stdio.h>
#include <errno.h>
#include <fitsio.h>

#include "splinetable.h"

static int parsefitstable(fitsfile *fits, struct splinetable *table);

int
readsplinefitstable(const char *path, struct splinetable *table)
{
	fitsfile *fits;
	int error = 0;

	fits_open_file(&fits,path,READONLY, &error);
	if (error != 0)
		return (error);

	error = parsefitstable(fits,table);
	fits_close_file(fits, &error);
	fits_report_error(stderr, error);

	return (error);
}

static int
parsefitstable(fitsfile *fits, struct splinetable *table)
{
	int error = 0;
	int hdus, type, i;
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

	fits_read_key(fits, TINT, "ORDER", &table->order, NULL, &error);
	table->periods = malloc(sizeof(table->periods[i])*table->ndim);
	for (i = 0; i < table->ndim; i++) {
		char name[255];
		sprintf(name,"PERIOD%d",i);
		fits_read_key(fits, TDOUBLE, name, &table->periods[i], NULL, &error);
	}

	/*
	 * Read the coefficient table
	 */

	table->naxes = malloc(sizeof(long)*table->ndim);
	fits_get_img_size(fits, table->ndim, table->naxes, &error);

	arraysize = table->naxes[0];
	for (i = 1; i < table->ndim; i++)
		arraysize *= table->naxes[i];
	table->coefficients = malloc(sizeof(double)*arraysize);

	fpixel = malloc(sizeof(long)*table->ndim);
	for (i = 0; i < table->ndim; i++)
		fpixel[i] = 1;

	fits_read_pix(fits, TDOUBLE, fpixel, arraysize, NULL,
	    table->coefficients, NULL, &error);

	free(fpixel);

	if (error != 0) {
		free(table->naxes);
		free(table->coefficients);
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
		if (error != 0)
			break;

		table->knots[i] = malloc(sizeof(double)*table->nknots[i]);
		fits_read_pix(fits, TDOUBLE, &fpix, table->nknots[i], NULL,
		    table->knots[i], NULL, &error);
	}

	return (error);
}

