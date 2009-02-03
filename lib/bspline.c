#include <math.h>

/*
 * Compute the value of the ith nth-order basis spline of a set
 * defined by knots at the point x.
 *
 * This is implemented using the De Boor algorithm, as outlined on
 * Wikipedia.
 */

double
bspline(double *knots, double x, int i, int n)
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
splineeval(double *knots, int nknots, double x, int order, int center)
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

	while (i < nknots-order-1 && i < center + order) {
		work += bspline(knots, x, i, order);
		i++;
	}

	return work;
}

