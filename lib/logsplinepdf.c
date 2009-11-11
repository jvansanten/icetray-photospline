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
	double mint, maxt;

	/* Find the boundaries by choosing the first and last points
	 * with full support */
	mint = table->knots[dim][table->order[dim]];
	maxt = table->knots[dim][table->nknots[dim] - table->order[dim] - 3];

	do {
		/*
		 * Get a starting point from the proposal distribution,
		 * making sure the PDF at the starting point is finite
		 * to avoid numerical problems.
		 */

		coords[dim] = lastval = (*proposal)();

		lastproppdf = (*proposal_pdf)(lastval,lastval);
		tablesearchcenters(table, coords, centers);
		lastlogpdf = ndsplineeval(table, coords, centers, 0);
		if (derivatives)
			lastlogpdf += log(ndsplineeval(table, coords, centers,
			    derivatives));
		accepted = 0;
	} while (!isfinite(lastlogpdf) || lastval < mint || lastval > maxt);

	for (i = -burnin; i < results; i++) {
		coords[dim] = val = (*proposal)();

		/*
		 * If we ended up outside the table, reject the sample
		 * and try again
		 */
		if (val > maxt || val < mint ||
		    tablesearchcenters(table, coords, centers) != 0) 
			goto reject;
			
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
		 * period has elapsed, write it to output array.
		 *
		 * NB: This code is used in both the accept and reject
		 * cases.
		 */

	     reject:
		if (i >= 0)
			result[i] = lastval;
	}
	#ifdef DEBUG
		printf("Efficiency: %e\n", (double)(accepted)/(double)(i));
	#endif
}

