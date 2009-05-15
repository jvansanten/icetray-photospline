#ifndef GLAM_H
#define GLAM_H

#include <suitesparse/cholmod.h>

#include "splinetable.h"

void glamfit(struct ndsparse *data, double *weights, double **coords,
    struct splinetable *out, double smooth, int order, int verbose,
    cholmod_common *c);

#endif
