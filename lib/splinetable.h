#ifndef _SPLINE_TABLE_H
#define _SPLINE_TABLE_H

struct splinetable {
	int ndim;
	int order;
	
	double **knots;
	long *nknots;

	double *periods;

	double *coefficients;
	long *naxes;
};

int readsplinefitstable(char *path, struct splinetable *table);

#endif /* _SPLINE_TABLE_H */

