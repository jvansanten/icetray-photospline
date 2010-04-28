#ifndef SPLINEUTIL_H
#define SPLINEUTIL_H

#include <suitesparse/cholmod.h>

struct ndsparse {
	/*
	 * This is an ntuple, and is similar to the CHOLMOD triplet
	 * type.
	 */

	size_t rows;
	size_t ndim;
	
	/* i contains coordinates, x the data values */
	double *x;
	int **i;
	int *ranges; /* Index ranges for each column */
};

cholmod_sparse *bsplinebasis(double *knots, size_t nknots, double *x,
    size_t npts, int order, cholmod_common *c);

int slicemultiply(struct ndsparse *a, cholmod_sparse *b, int dim,
    cholmod_common *c);
cholmod_sparse *kronecker_product(cholmod_sparse *a, cholmod_sparse *b,
    cholmod_common *c);
cholmod_sparse *box(cholmod_sparse *a, cholmod_sparse *b, cholmod_common *c);
cholmod_sparse *cholmod_tril(int dim, cholmod_common *c);

cholmod_dense *nnls_normal_block(cholmod_sparse *AtA, cholmod_dense *Atb,
   int verbose, cholmod_common *c);

void print_sparse(cholmod_sparse *a, cholmod_common *c);

#endif // SPLINEUTIL_H

