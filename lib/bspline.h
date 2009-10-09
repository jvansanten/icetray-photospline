#ifndef _BSPLINE_H
#define _BSPLINE_H

#include "splinetable.h"

/*
 * Compute the value of the ith nth-order basis spline of a set
 * defined by knots at the point x.
 */

double bspline(const double *knots, double x, int i, int n);
double bspline_deriv(const double *knots, double x, int i, int n);

/*
 * Evaluates the results of a full spline basis given a set of knots,
 * a position, an order, and a central spline for the position (or -1).
 * The central spline should be the index of the 0th order basis spline
 * that is non-zero at the position x.
 */

double splineeval(double *knots, double *weights, int nknots, double x,
    int order, int center);

/*
 * Spline table based hypersurface evaluation. ndsplineeval() takes a spline
 * coefficient table, a vector at which to evaluate the surface, and a vector
 * indicating the evaluation centers, as for splineeval().
 *
 * tablesearchcenters() provides a method to acquire a centers vector
 * for ndsplineeval() using a binary search. Depending on how the table
 * was produced, a more efficient method may be available.
 */

int tablesearchcenters(struct splinetable *table, double *x, int *centers);

double ndsplineeval(struct splinetable *table, const double *x, 
    const int *centers);


#endif /* _BSPLINE_H */
