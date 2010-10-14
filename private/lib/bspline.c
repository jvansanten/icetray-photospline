#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "photospline/bspline.h"

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

/*
 * A brain-dead reimplementation of de Boor's BSPLVB, which generates
 * the values of the non-zero B-splines at x from the bottom up without
 * unnccessarily recalculating terms. 
 * 
 * NB: the indexing used here assumes that x is fully supported, i.e. has
 * jhigh knots to the left and right. If this is not the case, the output
 * will be a casserole of nonsense.
 *
 * See Chapter X in: 
 * 
 * Carl de Boor. A Practical Guide to Splines, volume 27 of Applied
 *     Mathematical Sciences. Springer-Verlag, 1978.
 */

void
bsplvb_simple(const double *knots, double x, int left, int jhigh,
    float *restrict biatx)
{
	int i, j;
	double saved, term;
	double delta_l[jhigh], delta_r[jhigh];
	
	biatx[0] = 1.0;
		
	for (j = 0; j < jhigh-1; j++) {
		delta_r[j] = knots[left+j+1] - x;
		delta_l[j] = x - knots[left-j];
		
		saved = 0.0;
		
		for (i = 0; i < j+1; i++) {
			term = biatx[i] / (delta_r[i] + delta_l[j-i]);
			biatx[i] = saved + delta_r[i]*term;
			saved = delta_l[j-i]*term;
		}
		
		biatx[j+1] = saved;
	}
}

void
bsplvb(const double *knots, const double x, const int left, const int jlow,
    const int jhigh, float *restrict biatx,
    double *restrict delta_l, double *restrict delta_r)
{
	int i, j;
	double saved, term;

	if (jlow == 0)
		biatx[0] = 1.0;
		
	for (j = jlow; j < jhigh-1; j++) {
		delta_r[j] = knots[left+j+1] - x;
		delta_l[j] = x - knots[left-j];
		
		saved = 0.0;
		
		for (i = 0; i < j+1; i++) {
			term = biatx[i] / (delta_r[i] + delta_l[j-i]);
			biatx[i] = saved + delta_r[i]*term;
			saved = delta_l[j-i]*term;
		}
		
		biatx[j+1] = saved;
	}
}


void
bspline_deriv_nonzero(const double *knots, const double x, const int left, 
    const int n, float *restrict biatx)
{
	int i;
	double temp, a;
	
	/* Special case for constant splines */
	if (n == 0)
		return;
	
	/* Get the non-zero n-1th order B-splines at x */
	bsplvb_simple(knots, x, left, n, biatx);
	
	/* 
	 * Now, form the derivatives of the nth order B-splines from
	 * linear combinations of the lower-order splines.
	 */
	
	/* 
	 * On the last supported segment of the ith nth order spline,
	 * only the i+1th n-1th order spline is nonzero.
	 */
	temp = biatx[0];
	biatx[0] =  - n*temp / ((knots[left+1] - knots[left+1-n]));
	
	/* On the middle segments, both the ith and i+1th splines contribute. */
	for (i = 1; i < n; i++) {
		a = n*temp/((knots[left+i] - knots[left+i-n]));
		temp = biatx[i];
		biatx[i] = a - n*temp/(knots[left+i+1] - knots[left+i+1-n]);
	}
	/*
	 * On the first supported segment of the i+nth nth order spline,
	 * only the ith n-1th order spline is nonzero.
	 */
	biatx[n] = n*temp/((knots[left+n] - knots[left]));
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

	result = n * bspline(knots, x, i, n-1) / (knots[i+n] - knots[i]);
	result -= n * bspline(knots, x, i+1, n-1) / (knots[i+n+1] - knots[i+1]);
	
	return result;
}

double
bspline_deriv_2(const double *knots, double x, int i, int n)
{
	double result;

	if (n <= 1) {
		/*
		 * Special case the 1st order case, where B-Splines
		 * are linear functions from one knot to the next.
		 */

		return 0.0;
	}
	
	result = bspline(knots, x, i, n-2) /
	    ((knots[i+n] - knots[i])*(knots[i+n-1] - knots[i]));
	result -= bspline(knots, x, i+1, n-2) *
	    (1./(knots[i+n] - knots[i]) + 1./(knots[i+n+1] - knots[i+1])) / 
	    (knots[i+n] - knots[i+1]);
	result += bspline(knots, x, i+2, n-2) / 
	    ((knots[i+n+1] - knots[i+1])*(knots[i+n+1] - knots[i+2]));
	
	result *= n*(n-1);
	
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
	int i, min, max;

	for (i = 0; i < table->ndim; i++) {

		/*
		 * Do some sanity checks. Even inside the table, the results
		 * can make no sense (or worse, crash) if we are only
		 * a few knots in due to partial support.
		 */

		if (x[i] < table->knots[i][table->order[i]] ||
		    x[i] > table->knots[i][table->naxes[i]])
			return (-1);

		min = table->order[i];
		max = table->nknots[i]-2;
		do {
			centers[i] = (max+min)/2;

			if (x[i] < table->knots[i][centers[i]])
				max = centers[i]-1;
			else
				min = centers[i]+1;
		} while (x[i] < table->knots[i][centers[i]] ||
		    x[i] >= table->knots[i][centers[i]+1]);

		/*
		 * B-splines are defined on a half-open interval. For the
		 * last point of the interval, move center one point to the
		 * left to get the limit of the sum without evaluating
		 * absent basis functions.
		 */
		if (centers[i] == table->naxes[i])
			centers[i]--;
	}


	return (0);
}

static int
maxorder(int *order, int ndim)
{
	int i, max = 0;
	
	for (i = 0; i < ndim; i++)
		if (order[i] > max)
			max = order[i];
	
	return (max);
}
   
static float 
localbasis_sub(const struct splinetable *table, const int *centers,
    int n, int *restrict pos, unsigned long stride,
    const float *restrict localbasis[table->ndim])
{
	float acc = 0.0;
	int k;

	if (__builtin_expect(n+1 == table->ndim, 0)) {
		/*
		 * If we are at the last recursion level, the weights are
		 * linear in memory, so grab the row-level basis functions
		 * and multiply by the weights. Hopefully the compiler
		 * vector optimizations pick up on this code.
		 */

		const int woff = stride + centers[n] - table->order[n];
		const float *row = localbasis[n];
		for (k = 0; k <= table->order[n]; k++)
			acc += table->coefficients[woff+k]*row[k];
	} else {
		for (k = -table->order[n]; k <= 0; k++) {
			/*
			 * If we are not at the last dimension, record where we
			 * are, multiply in the row basis value, and recurse.
			 */

			pos[n] = centers[n] + k;
			acc += localbasis_sub(table, centers, n+1, pos,
			    stride + pos[n]*table->strides[n], localbasis) *
			    localbasis[n][k+table->order[n]];
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
ndsplineeval(struct splinetable *table, const double *x, const int *centers,
    int derivatives)
{
	int n, offset;
	int maxdegree = maxorder(table->order, table->ndim) + 1; 
	int pos[table->ndim];
	double result;
	float localbasis[table->ndim][maxdegree];
	const float *localbasis_ptr[table->ndim];
	
	for (n = 0; n < table->ndim; n++) {
		if (derivatives & (1 << n)) {
			bspline_deriv_nonzero(table->knots[n], x[n], centers[n],
			    table->order[n], localbasis[n]);
		} else {
			bsplvb_simple(table->knots[n], x[n], centers[n],
			    table->order[n] + 1, localbasis[n]);
		}

		localbasis_ptr[n] = localbasis[n];
	}

	result = localbasis_sub(table, centers, 0, pos, 0, localbasis_ptr);

	return (result);
}

double
ndsplineeval_deriv2(struct splinetable *table, const double *x, const int *centers,
    int derivatives)
{
	int n, offset;
	int maxdegree = maxorder(table->order, table->ndim) + 1; 
	int pos[table->ndim];
	double result;
	float localbasis[table->ndim][maxdegree];
	const float *localbasis_ptr[table->ndim];
	
	for (n = 0; n < table->ndim; n++) {
		if (derivatives & (1 << n)) {
			for (offset = -table->order[n]; offset <= 0; offset++) {
				localbasis[n][offset+table->order[n]] =
				    bspline_deriv_2(table->knots[n],x[n],
				    centers[n] + offset, table->order[n]);
			}
		} else {
			bsplvb_simple(table->knots[n], x[n], centers[n],
			    table->order[n] + 1, localbasis[n]);
		}

		localbasis_ptr[n] = localbasis[n];
	}

	result = localbasis_sub(table, centers, 0, pos, 0, localbasis_ptr);

	return (result);
}


