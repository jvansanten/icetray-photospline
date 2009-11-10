#include <math.h>
#include <gsl/gsl_rng.h>

#include "splinetable.h"
#include "bspline.h"
#include "logsplinepdf.h"

void logsplinepdf_n_sample(double *result, int results, int burnin,
    double *coords, int dim, struct splinetable *table, int derivatives,
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
	lastlogpdf = ndsplineeval(table, coords, centers, 0);
	if (derivatives)
		lastlogpdf += log(ndsplineeval(table, coords, centers,
		    derivatives));
	accepted = 0;

	for (i = -burnin; i < results; i++) {
		coords[dim] = val = (*proposal)();
		if (tablesearchcenters(table, coords, centers) != 0) {
			// If we ended up outside the table, try again
			i--; continue;
		}
			
		logpdf = ndsplineeval(table, coords, centers, 0);
		if (derivatives)
			logpdf += log(ndsplineeval(table, coords, centers,
			    derivatives));
		proppdf = (*proposal_pdf)(val,lastval);
		odds = exp(logpdf - lastlogpdf);
		odds *= lastproppdf/proppdf;

		if (odds > 1. || gsl_rng_uniform(rng) < odds) {
			/* Accept this value */
			lastval = val;
			lastlogpdf = logpdf;
			lastproppdf = proppdf;
			#ifdef DEBUG
				accepted++;
			#endif
		}

		/*
		 * Lastval has whatever we decided on now. If the burn-in
		 * period has elapsed, write it to output array
		 */

		if (i >= 0)
			result[i] = lastval;
	}
	#ifdef DEBUG
		printf("Efficiency: %e\n", (double)(accepted)/(double)(i));
	#endif
}

