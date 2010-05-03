#include <stdio.h>
#include <string.h>
#include <math.h>

#include <time.h>

#include <suitesparse/cholmod.h>
#include "splineutil.h"

#define MAX_TRIALS 10

static int intcmp(const void *xa, const void *xb);

/*
 * An implementation of the Portugal/Judice/Vicente block-pivoting algorithm for
 * pre-formulated normal equations
 *
 *    See:
 *       A Comparison of Block Pivoting and Interior-Point Algorithms for
 *       Linear Least Squares Problems with Nonnegative Variables
 *
 *       Author(s): Luis F. Portugal, Joaquim J. Judice, Luis N. Vicente
 *       Source: Mathematics of Computation, Vol. 63, No. 208 (Oct., 1994),
 *         pp. 625-643
 *       Published by: American Mathematical Society
 *       Stable URL: http://www.jstor.org/stable/2153286
 */

cholmod_dense *
nnls_normal_block(cholmod_sparse *AtA, cholmod_dense *Atb, int verbose,
   cholmod_common *c)
{
	int nvar = AtA->nrow;
	long F[nvar], G[nvar];
	int H1[nvar], H2[nvar];
	cholmod_dense *x, *y, *x_F, *Atb_F;
	cholmod_sparse *AtA_F;
	cholmod_factor *L;
	int nF, nG, nH1, nH2, ninf;
	int i, j, k, trials, iter;

	clock_t t0,t1;

	/* XXX: make these settable? */
	iter = 3*nvar;		/* Maximum number of iterations */
	trials = MAX_TRIALS;	/* Runs without progress before reverting
				 * to a deterministic algorithm */

	x = cholmod_l_zeros(nvar, 1, CHOLMOD_REAL, c);
	y = cholmod_l_allocate_dense(nvar, 1, nvar, CHOLMOD_REAL, c);

	/*
	 * Initialize variables such that all coefficients are in the
	 * active set.
	 */

	nF = 0;			/* No feasible coefficients */
	nG = nvar;
	for (i = 0; i < nvar; i++)
		G[i] = i;
	ninf = nvar + 1;	/* Number of infeasible coefficients */

	/* Set up the dual vector */
	for (i = 0; i < nvar; i++)
		((double *)(y->x))[i] = -((double *)(Atb->x))[i];

	while (iter-- > 0) {
		/*
		 * Fill H1, H2 with the points that must be moved from the
		 * passive to active and active to passive sets, respectively.
		 */

		nH1 = nH2 = 0;
		for (i = 0; i < nF; i++)
			if (((double *)(x->x))[F[i]] < 0)
				H1[nH1++] = F[i];
		for (i = 0; i < nG; i++)
			if (((double *)(y->x))[G[i]] < 0)
				H2[nH2++] = G[i];

		/*
		 * If all coefficients were positive, then we are done.
		 */

		if (nH1 == 0 && nH2 == 0)
			break;

		/*
		 * Check the status of the bulk set switching.
		 */
		
		if (nH2 + nH1 < ninf) {
			ninf = nH2 + nH1;
			trials = MAX_TRIALS;
		} else {
			/* Stuck, check if we need to try something else */
			trials--;
			if (trials < 0) {
				/*
				 * Out of luck. Fall back to slow but
				 * guaranteed method (picking the last
				 * infeasible coordinate).
				 */
				
				if (H1[nH1 - 1] > H2[nH1 - 2]) {
					H1[0] = H1[nH1 - 1];
					nH1 = 1; nH2 = 0;
				} else {
					H2[0] = H2[nH2 - 1];
					nH2 = 1; nH1 = 0;
				}
			}
		}

		if (verbose)
			printf("Infeasibles: %d\n", ninf);

		/*
		 * Next, remove elements in H1 from F, and add them to G,
		 * exploiting the fact that H1 elements are in order
		 * relative to their order in F.
		 */
		for (i = 0, j = 0; i < nH1; i++) {
			G[nG++] = H1[i];
			while (F[j] != H1[i]) j++;
			for (k = j+i; k+1 < nF; k++)
				F[k-i] = F[k-i+1];
		}
		nF -= nH1;

		/* And vice versa */
		for (i = 0, j = 0; i < nH2; i++) {
			F[nF++] = H2[i];
			while (G[j] != H2[i]) j++;
			for (k = j+i; k+1 < nG; k++)
				G[k-i] = G[k-i+1];
		}
		nG -= nH2;

		qsort(G, nG, sizeof(G[0]), intcmp);
		qsort(F, nF, sizeof(F[0]), intcmp);

		/* Solve the unconstrained part */

		if (verbose)
		    printf("Unconstrained solve for %d of %d coefficients\n",
		      nF, nvar);
		AtA_F = cholmod_l_submatrix(AtA, F, nF, F, nF, 1, 1, c);
		AtA_F->stype = 1;
		Atb_F = cholmod_l_allocate_dense(nF, 1, nF, CHOLMOD_REAL, c);
		for (i = 0; i < nF; i++)
			((double *)(Atb_F->x))[i] = ((double *)(Atb->x))[F[i]];

		t0 = clock();

		L = cholmod_l_analyze(AtA_F, c);

		t1 = clock();
		printf("Analyze: %f s\n",(double)(t1-t0)/(CLOCKS_PER_SEC));
		t0 = clock();

		cholmod_l_factorize(AtA_F, L, c);

		t1 = clock();
		printf("Factorize: %f s\n",(double)(t1-t0)/(CLOCKS_PER_SEC));
		t0 = clock();

		x_F = cholmod_l_solve(CHOLMOD_A, L, Atb_F, c);

		t1 = clock();
		printf("Solve: %f s\n",(double)(t1-t0)/(CLOCKS_PER_SEC));
		t0 = clock();

		cholmod_l_free_factor(&L, c);
		for (i = 0; i < nF; i++)
			((double *)(x->x))[F[i]] = ((double *)(x_F->x))[i];
		cholmod_l_free_sparse(&AtA_F, c);

		/* Update the constrained part */

		for (i = 0; i < nG; i++)
			((double *)(x->x))[G[i]] = 0;
		for (i = 0; i < nF; i++)
			((double *)(y->x))[F[i]] = 0;

		{
			cholmod_sparse *AtA_FG;
			cholmod_dense *Atb_G;
			double ones[2] = {1., 0}, mones[2] = {-1., 0};

			Atb_G = cholmod_l_allocate_dense(nG, 1, nG,
			    CHOLMOD_REAL, c);
			for (i = 0; i < nG; i++)
				((double *)(Atb_G->x))[i] =
				    ((double *)(Atb->x))[G[i]];

			AtA_FG = cholmod_l_submatrix(AtA, G, nG, F, nF, 1,1,c);
			cholmod_l_sdmult(AtA_FG, 0, ones, mones, x_F, Atb_G, c);
			
			for (i = 0; i < nG; i++)
				((double *)(y->x))[G[i]] =
				    ((double *)(Atb_G->x))[i];

			cholmod_l_free_dense(&Atb_G, c);
			cholmod_l_free_sparse(&AtA_FG, c);
		}
		
		cholmod_l_free_dense(&x_F, c);
	}

	cholmod_l_free_dense(&y, c);

	return (x);
}

static int
intcmp(const void *xa, const void *xb)
{
	const int *a, *b;
	a = xa; b = xb;

	if (*a < *b)
		return (-1);
	else if (*a > *b)
		return (1);
	else
		return (0);
}

