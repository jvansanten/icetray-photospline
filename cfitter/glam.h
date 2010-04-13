#ifndef GLAM_H
#define GLAM_H

#include <suitesparse/cholmod.h>

#include "splinetable.h"
#include "splineutil.h"

void glamfit(struct ndsparse *data, double *weights, double **coords,
    struct splinetable *out, double smooth, int *order, int *penorder,
    int verbose, cholmod_common *c);
    
void
glamfit_complex(struct ndsparse *data, double *weights, double **coords,
    struct splinetable *out, int *order, cholmod_sparse* penalty,
    int verbose, cholmod_common *c);

cholmod_sparse* add_penalty_term(long *nsplines, double *knots, int ndim, int dim, int order,
    int porder, double scale, 
    cholmod_sparse *penalty, cholmod_common *c);

#endif
