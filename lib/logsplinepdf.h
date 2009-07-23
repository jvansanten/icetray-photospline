#ifndef LOGSPLINEPDF
#define LOGSPLINEPDF

#include <gsl/gsl_rng.h>

#include "splinetable.h"

void logsplinepdf_n_sample(double *result, int results, int burnin,
    double *coords, int dim, struct splinetable *table,
    double (* proposal)(void), double (* proposal_pdf)(double, double),
    gsl_rng *rng);

#endif
