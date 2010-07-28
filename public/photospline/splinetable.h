#ifndef _SPLINE_TABLE_H
#define _SPLINE_TABLE_H

struct splinetable {
	int ndim;
	int *order;

	double **knots;
	long *nknots;

	double **extents;

	double *periods;

	double *coefficients;
	long *naxes;
	unsigned long *strides;

	int naux;
	char ***aux;
};

typedef enum {
	SPLINETABLE_INT,
	SPLINETABLE_DOUBLE
} splinetable_dtype;

int readsplinefitstable(const char *path, struct splinetable *table);
void splinetable_free(struct splinetable *table);
char * splinetable_get_key(struct splinetable *table, const char *key);
int splinetable_read_key(struct splinetable *table, splinetable_dtype type,
    char *key, void *result);


#endif /* _SPLINE_TABLE_H */

