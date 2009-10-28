#ifndef GLAM_H
#define GLAM_H

#include <suitesparse/cholmod.h>

#include "splinetable.h"
#include "splineutil.h"

void glamfit(struct ndsparse *data, double *weights, double **coords,
    struct splinetable *out, double smooth, int *order, int *penorder,
    int verbose, cholmod_common *c);

#endif
