#ifndef _BSPLINE_H
#define _BSPLINE_H

/*
 * Compute the value of the ith nth-order basis spline of a set
 * defined by knots at the point x.
 */

double bspline(const double *knots, double x, int i, int n);

/*
 * Evaluates the results of a full spline basis given a set of knots,
 * a position, an order, and a central spline for the position (or -1).
 * The central spline should be the index of the 0th order basis spline
 * that is non-zero at the position x.
 */

double splineeval(double *knots, double *weights, int nknots, double x, int order,
    int center);

double ndsplineeval(double **knots, const double *weights, int ndim, 
    const int *nknots, const double *x, int order, const int *centers);


#endif /* _BSPLINE_H */
