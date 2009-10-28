#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bspline.h"

/*
 * FASTANDLOOSE mode avoids some expensive malloc() calls
 * at the price of loss of generality (tables with more than
 * MAXORDER dimensions cannot be processed in this mode).
 */

#define	FASTANDLOOSE	1
#define	MAXORDER	5

/*
 * Compute the value of the ith nth-order basis spline of a set
 * defined by knots at the point x.
 *
 * This is implemented using the De Boor algorithm, as outlined on
 * Wikipedia.
 */

double
bspline(const double *knots, double x, int i, int n)
{
	double result;

	if (n == 0) {
		/*
		 * Special case the 0th order case, where B-Splines
		 * are constant functions from one knot to the next.
		 */

		if (x >= knots[i] && x < knots[i+1])
			return 1.0;
		else
			return 0.0;
	}

	result = (x - knots[i])*bspline(knots, x, i, n-1) /
	    (knots[i+n] - knots[i]);
	result += (knots[i+n+1] - x)*bspline(knots, x, i+1, n-1) /
	    (knots[i+n+1] - knots[i+1]);

	return result;
}

double
bspline_deriv(const double *knots, double x, int i, int n)
{
	double result;

	if (n == 0) {
		/*
		 * Special case the 0th order case, where B-Splines
		 * are constant functions from one knot to the next.
		 */

		return 0.0;
	}

	result = ((x - knots[i])*bspline_deriv(knots, x, i, n-1) + 
	    bspline(knots, x, i, n-1)) / (knots[i+n] - knots[i]);
	result += ((knots[i+n+1] - x)*bspline_deriv(knots, x, i+1, n-1) -
	    bspline(knots, x, i+1, n-1)) / (knots[i+n+1] - knots[i+1]);

	return result;
}

/*
 * Evaluates the results of a full spline basis given a set of knots,
 * a position, an order, and a central spline for the position (or -1).
 * The central spline should be the index of the 0th order basis spline
 * that is non-zero at the position x.
 */

double
splineeval(double *knots, double *weights, int nknots, double x, int order,
    int center)
{
	double work = 0.0;
	int i;

	if (center < 0) {
		/* XXX: should be a binary search */
		for (center = 0; center+1 < nknots; center++) {
			if (x > knots[center] && x < knots[center+1])
				break;
		}
	
		if (center+1 >= nknots)
			return 0.0;
	}

	i = center - order;
	if (i < 0)
		i = 0;

	while (i < nknots-order-1 && i <= center) {
		work += weights[i]*bspline(knots, x, i, order);
		i++;
	}

	return work;
}

int
tablesearchcenters(struct splinetable *table, double *x, int *centers)
{
	int i;

	for (i = 0; i < table->ndim; i++) {

		/*
		 * Do some sanity checks. Even inside the table, the results
		 * can make no sense (or worse, crash) if we are only
		 * a few knots in due to partial support.
		 */

		if (x[i] < table->knots[i][table->order[i]] ||
		    x[i] > table->knots[i][table->nknots[i]-table->order[i]-1])
			return (-1);

		/* XXX: should be a binary search */
		for (centers[i] = table->order[i];
		    centers[i]+1 < table->nknots[i]; centers[i]++) {
			if (x[i] >= table->knots[i][centers[i]] &&
			    x[i] < table->knots[i][centers[i]+1])
				break;
		}
		if (centers[i]+2 >= table->nknots[i])
			return (-1);
	}


	return (0);
}
   
static double
localbasis_sub(const double *weights, const int *centers, int ndim,
    int *order, int n, const long *naxes, const unsigned long *strides,
#if FASTANDLOOSE
    int pos[ndim], unsigned long stride, double localbasis[ndim][MAXORDER])
#else
    int pos[ndim], unsigned long stride, double *localbasis[ndim])
#endif
{
	double acc = 0.0;
	int k;

	if (n+1 == ndim) {
		/*
		 * If we are at the last recursion level, the weights are
		 * linear in memory, so grab the row-level basis functions
		 * and multiply by the weights. Hopefully the compiler
		 * vector optimizations pick up on this code.
		 */

		long woff;
		woff = stride + centers[n] - order[n];

		for (k = 0; k <= order[n]; k++) {
			acc += weights[k + woff]*
			    localbasis[n][k];
		}
	} else {
		for (k = -order[n]; k <= 0; k++) {
			/*
			 * If we are not at the last dimension, record where we
			 * are, multiply in the row basis value, and recurse.
			 */

			pos[n] = centers[n] + k;
			acc += localbasis_sub(weights, centers, ndim, order,
			    n+1, naxes, strides, pos,
			    stride + pos[n]*strides[n], localbasis)
			    * localbasis[n][k+order[n]];
		}
	}

	return acc;
}

/*
 * The N-Dimensional tensor product basis version of splineeval.
 * Evaluates the results of a full spline basis given a set of knots,
 * a position, an order, and a central spline for the position (or -1).
 * The central spline should be the index of the 0th order basis spline
 * that is non-zero at the position x.
 *
 * x is the vector at which we will evaluate the space
 */

double
ndsplineeval(struct splinetable *table, const double *x, const int *centers)
{
	int n, offset;
	int pos[table->ndim];
	double result;
	#if FASTANDLOOSE
	double localbasis[table->ndim][MAXORDER];
	#else
	double *localbasis[table->ndim];
	#endif

	for (n = 0; n < table->ndim; n++) {
		#if FASTANDLOOSE
		   if (table->order[n] > MAXORDER) {
			fprintf(stderr, "B-spline FASTANDLOOSE mode is on "
			    "and the table order along dimension %d (%d) is "
			    "larger than MAXORDER (%d). Either disable FASTAND"
			    "LOOSE mode or adjust MAXORDER in bspline.c.\n",
			    n, table->order[n], MAXORDER);
			abort();
		   }
		#else
		   localbasis[n] = calloc(table->order[n] + 1, sizeof(double));
		#endif
		for (offset = -table->order[n]; offset <= 0; offset++) {
			localbasis[n][offset+table->order[n]] =
			     bspline(table->knots[n],x[n],
			         centers[n] + offset, table->order[n]);
		}
	}

	result = localbasis_sub(table->coefficients, centers, table->ndim,
	    table->order, 0, table->naxes, table->strides, pos, 0, localbasis);

	#if !(FASTANDLOOSE)
	for (n = 0; n < table->ndim; n++)
		free(localbasis[n]);
	#endif

	return (result);
}
