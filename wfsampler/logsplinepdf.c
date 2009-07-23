#include <math.h>
#include <gsl/gsl_rng.h>

#include "splinetable.h"
#include "bspline.h"
#include "logsplinepdf.h"

void logsplinepdf_n_sample(double *result, int results, int burnin,
    double *coords, int dim, struct splinetable *table,
    double (* proposal)(void), double (* proposal_pdf)(double, double),
    gsl_rng *rng)
{
	int i, accepted;
	int centers[table->ndim];
	double val, lastval;
	double lastlogpdf, logpdf;
	double lastproppdf, proppdf;
	double odds, test;

	coords[dim] = lastval = (*proposal)();
	lastproppdf = (*proposal_pdf)(lastval,lastval);
	tablesearchcenters(table, coords, centers);
	lastlogpdf = ndsplineeval(table, coords, centers);

	for (i = -burnin; i < results; i++) {
		coords[dim] = val = (*proposal)();
		tablesearchcenters(table, coords, centers);
		logpdf = ndsplineeval(table, coords, centers);
		proppdf = (*proposal_pdf)(val,lastval);
		odds = exp(logpdf - lastlogpdf);
		odds *= lastproppdf/proppdf;

		if (odds > 1. || gsl_rng_uniform(rng) < odds) {
			/* Accept this value */
			lastval = val;
			lastlogpdf = logpdf;
			lastproppdf = proppdf;
			if (i >= 0)
				result[i] = val;
			accepted++;
		} else {
			if (i >= 0)
				result[i] = lastval;
		}
	}
	
	printf("Efficiency: %e\n", (double)(accepted)/(double)(results));
}

