#include <stdlib.h>
#include <float.h>
#include <string.h>

#include "photospline/bspline.h"

static double divdiff(double *x, double *y, size_t n)
{
	if (n == 1)
		return y[0];
	
	return ((divdiff(&x[1], &y[1], n-1) - divdiff(x, y, n))
	    / (x[n-1] - x[0]));
}

static int factorial(int n)
{
	int i = n-1;
	int acc = n;
	
	for ( ; i > 1; i--)
		acc *= i;
	
	return (acc);
}

static int
double_cmp(const void *xa, const void *xb)
{
	const double *a, *b;
	a = xa; b = xb;

	if (*a < *b)
		return (-1);
	else if (*a > *b)
		return (1);
	else
		return (0);
}

/* 
 * The local blossom of the convolution of the splines defined on knot
 * vectors x and y can be evaluated at point z via iterated divided differences.
 * 
 * This is analogous to Stroem Equation 13 and Lemma 9, but with the prefactor
 * adapted to account for the fact that one of the splines is de-Boor normalized
 * and the other unit normalized.
 */ 
double convoluted_blossom(double *x, size_t nx, double *y, size_t ny, double z,
    double *bags, size_t nbags)
{
	double scale, fun_x[nx], fun_y[ny];
	int i, j, k;
	
	scale = x[nx-1] - x[0];
	
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (x[i] + y[j] - z > 0) {
				double det = 1;
				for (k = 0; k < nbags; k++) {
					det *= (x[i] + y[j] - bags[k]);
				}
				fun_y[j] = det;
			} else {
				fun_y[j] = 0;
			}
		}
		fun_x[i] = divdiff(y, fun_y, ny);
	}
	return (scale*divdiff(x, fun_x, nx));
}

/* The coefficient of a spline in the convolved basis is a linear
 * combination of the blossoms of each spline in the un-convolved basis convolved
 * with the kernel spline. Here we just store the raw blossoms and multiply by
 * the coefficients of each un-convolved spline later.
 * 
 * This is analogous Stroem, Proposition 10, but with the prefactor adapted to
 * account for the fact that one of the splines is de-Boor normalized and the
 * other unit normalized
 */

int splinetable_convolve(struct splinetable *table, int dim, double *knots, size_t n_knots)
{
	double *rho, *coefficients, *bundle, **trafo, norm;
	size_t n_rho, arraysize;
	unsigned long *strides;
	long *naxes;
	unsigned convorder;
	int i, j, k, q;
	
	/* Construct the new knot field. */
	rho = calloc(table->nknots[dim], sizeof(double));
	n_rho = 0;
	for (i = 0; i < table->nknots[dim]; i++)
		for (j = 0; j < n_knots; j++)
			rho[n_rho++] = table->knots[dim][i] + knots[j];
	
	/* Order the new knot field and remove any duplicates. */
	qsort(rho, n_rho, sizeof(rho[0]), double_cmp);
	for (i = 1; i < n_rho; i++) {
		if (rho[i] - rho[i-1] < DBL_EPSILON) {
			for (j = i; j < n_rho; j++)
				rho[j] = rho[j+1];
			n_rho--;
		}
	}
	
	convorder = table->order[dim] + n_knots - 1;
	
	/* Allocate space for the convolved coefficients */
	naxes = malloc(sizeof(long)*table->ndim);
	strides = malloc(sizeof(unsigned long)*table->ndim);
	
	memcpy(naxes, table->naxes, sizeof(long)*table->ndim);
	naxes[dim] = n_rho - convorder - 1;
	
	strides[table->ndim - 1] = arraysize = 1;
	for (i = table->ndim-1; i >= 0; i--) {
		arraysize *= naxes[i];
		if (i > 0)
			strides[i-1] = arraysize;
	}
			
	coefficients = malloc(sizeof(double)*arraysize);

	k = table->order[dim] + 1;
	q = convorder + 1;
	/*
	 * NB: we're convolving a de-Boor spline with a unit-norm spline,
	 * hence q!(k-1)! rather than (q-1)!(k-1)! (as for two de-Boor splines).
	 */
	norm = ((double)(factorial(q)*factorial(k-1)))/((double)factorial(k+q-1));
	if (k % 2 != 0)
		norm *= -1;
	
	/*
	 * Now, calculate a transformation from coefficients on the raw knot field
	 * to coefficients on the convoluted knot field. Since the knots are on a
	 * grid, this transformation can be applied to each slice of the array.
	 */
	trafo = malloc(sizeof(double*)*naxes[dim]);
	for (i = 0; i < naxes[dim]; i++) {
		trafo[i] = malloc(sizeof(double)*table->naxes[dim]);
		for (j = 0; j < table->naxes[dim]; j++) {
			trafo[i][j] = norm*convoluted_blossom(&table->knots[dim][j],
			    k+1, knots, n_knots, rho[i], &rho[i+1], k+q-1);
		}
	}
	
	/*
	 * Walk through slices along the chosen dimension, applying the
	 * transformation to each.
	 */
	size_t n_slices = arraysize/naxes[dim];
	size_t offset;
	for (i = 0; i < n_slices; i++) {
		for (j = 0; j < naxes[dim]; j++) {
			/* XXX TODO: figure out how to set up strides */
		}
	}
	
	/* Swap out the new components of the table */
	free(table->coefficients);
	free(table->naxes);
	free(table->strides);
	free(table->knots[dim]);
	
	table->coefficients = coefficients;
	table->naxes = naxes;
	table->strides = strides;
	table->knots[dim] = rho;
	
	table->nknots[dim] = n_rho;
	table->order[dim] = convorder;
	
	table->extents[dim][0] = table->knots[dim][convorder];
	table->extents[dim][1] = table->knots[dim][table->nknots[dim]-convorder-1];
	
	return (0);
}
