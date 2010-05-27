#include <stdio.h>
//#include <string.h>
#include <math.h>

#include <time.h>
#include <float.h>

#include <suitesparse/cholmod.h>
#include "splineutil.h"

cholmod_dense*
cholesky_solve(cholmod_sparse *AtA, cholmod_dense *Atb, cholmod_common *c, int verbose, int n_resolves)
{
	int i,j,nvar;
	double ones[2] = {1., 0}, mones[2] = {-1., 0};
	double sum_delta_b;
	cholmod_dense *x, *delta_x, *delta_Atb;
	cholmod_factor *L;
	clock_t t0, t1;

	/* drop small entries from AtA */
	cholmod_l_drop(DBL_EPSILON, AtA, c);

	/* set symmetry */
	AtA->stype = 1;

	/* allocate workspace */
	nvar = AtA->nrow;
	if (n_resolves > 0) {
		delta_Atb = cholmod_l_allocate_dense(nvar, 1, nvar,
                            CHOLMOD_REAL, c);
		delta_x = cholmod_l_allocate_dense(nvar, 1, nvar,
                            CHOLMOD_REAL, c);
	}
	
	t0 = clock();

	/* compute the nonzero pattern of the Cholesky factorization */	
	L = cholmod_l_analyze(AtA, c);

	if (verbose) {
		t1 = clock();
		printf("Analyze[%d]: %f s\n",nvar,(double)(t1-t0)/(CLOCKS_PER_SEC));
		t0 = clock();
	}

	/* compute the numerical values of the Cholesky factorization */	
	cholmod_l_factorize(AtA, L, c);

	if (verbose) {
		t1 = clock();
		printf("Factorize[%d]: %f s\n",nvar,(double)(t1-t0)/(CLOCKS_PER_SEC));
		t0 = clock();
	}

	/* solve the system AtA*x = Atb */
	x = cholmod_l_solve(CHOLMOD_A, L, Atb, c);

	if (verbose) {
		t1 = clock();
		printf("Solve[%d]: %f s\n",nvar,(double)(t1-t0)/(CLOCKS_PER_SEC));
		t0 = clock();
	}

	/* refine the solution if needed */
	for (j = 0; j < n_resolves; j++) {
		sum_delta_b = 0.;

		if (verbose) t0 = clock();

		/* copy the right-hand side */
		for (i = 0; i < nvar; i++)
		    ((double *)(delta_Atb->x))[i] =
		        ((double *)(Atb->x))[i];

		/* delta_b = A*(x + delta_x) - b */
		cholmod_l_sdmult(AtA, 0, ones, mones, x, delta_Atb, c);

		/* solve A*delta_x = delta_b */
		delta_x = cholmod_l_solve(CHOLMOD_A, L, delta_Atb, c);

		/* x' = x - delta_x */
		for (i = 0; i < nvar; i++) {
			sum_delta_b += (((double *)(delta_Atb->x))[i])*(((double *)(delta_Atb->x))[i]);
			((double *)(x->x))[i] -= ((double *)(delta_x->x))[i];
		}

		if (verbose) {
			t1 = clock();
			printf("reSolve %d: %f s (db: %e)\n",
		    	    j,(double)(t1-t0)/(CLOCKS_PER_SEC),sqrt(sum_delta_b));
		}

		cholmod_l_free_dense(&delta_x, c);
	}

	cholmod_free_factor(&L, c);
	if (n_resolves > 0) cholmod_free_dense(&delta_Atb, c);

	return (x);
}

