#include <math.h>
#include <stdio.h>

#include "bspline.h"

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

		if (x > knots[i] && x < knots[i+1])
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
		if (x[i] < table->knots[i][0])
			return (-1);

		/* XXX: should be a binary search */
		for (centers[i] = 0; centers[i]+1 < table->nknots[i];
		    centers[i]++) {
			if (x[i] > table->knots[i][centers[i]] &&
			    x[i] < table->knots[i][centers[i]+1])
				break;
		}

		if (centers[i]+1 >= table->nknots[i])
			return (-1);
	}

	return (0);
}
   
static double
localbasis_sub(const double *weights, const int *centers, int ndim,
    int order, int n, const long *naxes, int pos[ndim],
    double localbasis[ndim][order + 1])
{
	double acc = 0.0;
	int k;

	if (n+1 == ndim) {
		int i, j;
		int stride;

		/*
		 * If we are at the last recursion level, the weights are
		 * linear in memory, so grab the row-level basis functions
		 * and multiply by the weights. Hopefully the compiler
		 * vector optimizations pick up on this code.
		 */

		j = 0;
		stride = 1;
		for (i = n; i >= 0; i--) {
			if (i < n) j += pos[i]*stride;
			stride *= naxes[i];
		}

		for (k = -order; k <= 0; k++) {
			acc += weights[j + k + centers[n]]*
			    localbasis[n][k+order];
		}
	} else {
		for (k = -order; k <= 0; k++) {
			/*
			 * If we are not at the last dimension, record where we
			 * are, multiply in the row basis value, and recurse.
			 */

			pos[n] = centers[n] + k;
			acc += localbasis_sub(weights, centers, ndim, order,
			    n+1, naxes, pos, localbasis)*localbasis[n][k+order];
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
	double localbasis[table->ndim][table->order + 1];
	int pos[table->ndim];

	for (n = 0; n < table->ndim; n++) {
		for (offset = -table->order; offset <= 0; offset++) {
			localbasis[n][offset+table->order] =
			     bspline(table->knots[n],x[n],
			     centers[n] + offset, table->order);
		}
	}

	return (localbasis_sub(table->coefficients, centers, table->ndim,
	    table->order, 0, table->naxes, pos, localbasis));
}
