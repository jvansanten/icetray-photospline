#include <stdio.h>
#include <string.h>
#include <math.h>

#include <time.h>
#include <float.h>

#include <suitesparse/cholmod.h>
#include "splineutil.h"
#include "cholesky_solve.h"

#define MAX_TRIALS 3
#define N_RESOLVES 0
#define KKT_TOL 1e-6

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
	int nF, nG, nH1, nH2, ninf;
	int i, j, k, trials, murty_steps, iter;

	/* XXX: make these settable? */
	iter = 3*nvar;		/* Maximum number of iterations */
	trials = MAX_TRIALS;	/* Runs without progress before reverting
				 * to a deterministic algorithm */
	murty_steps = MAX_TRIALS;

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

	/* Drop small entries from the problem */
	//cholmod_l_drop(DBL_EPSILON, AtA, c);

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
			if (((double *)(x->x))[F[i]] < -KKT_TOL)
				H1[nH1++] = F[i];
		for (i = 0; i < nG; i++)
			if (((double *)(y->x))[G[i]] < -KKT_TOL)
				H2[nH2++] = G[i];

		/*
		 * If all coefficients were positive, then we are done.
		 */

		if (nH1 == 0 && nH2 == 0)
			break;

		if (ninf <= murty_steps)
			trials = -1;

		/*
		 * Check the status of the bulk set switching.
		 * After MAX_TRIALS iterations of Murty's finite method,
		 * revert to block switching.
		 */
		
		if (ninf > murty_steps &&
		    (nH2 + nH1 < ninf || trials < -murty_steps)) {
			if (nH2 + nH1 <= ninf)
			    murty_steps++;
			ninf = nH2 + nH1;
			trials = MAX_TRIALS;
		} else {
			/* Stuck, check if we need to try something else */
			trials--;
			if (verbose)
				printf("Stuck! trials: %d nH1: %d nH2: %d\n",trials,nH1,nH2);
			if (trials < 0) {
				/*
				 * Out of luck. Fall back to slow but
				 * guaranteed method (picking the last
				 * infeasible coordinate).
				 */
				
				if (nH2 == 0) {
					goto maxh1;
				} else if (nH1 == 0) {
					goto maxh2;
				} else if (H1[nH1 - 1] > H2[nH2 - 1]) {
					maxh1:
					H1[0] = H1[nH1 - 1];
					nH1 = 1; nH2 = 0;
					if (verbose)
						printf("H1: %d (%e)\n",H1[0],((double *)(x->x))[H1[0]]);
				} else {
					maxh2:
					H2[0] = H2[nH2 - 1];
					nH2 = 1; nH1 = 0;
					if (verbose)
						printf("H2: %d (%e)\n",H2[0],((double *)(y->x))[H2[0]]);
				}
			}
		}

		if (verbose)
			printf("Iteration %d Infeasibles: %d\n", (3*nvar - iter), ninf);

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

		/* Solve the system AtA_F*x = Atb_F, refining the solution iteratively */
		x_F = cholesky_solve(AtA_F, Atb_F, c, verbose, N_RESOLVES);

		for (i = 0; i < nF; i++)
			((double *)(x->x))[F[i]] = ((double *)(x_F->x))[i];
		cholmod_l_free_sparse(&AtA_F, c);
		cholmod_l_free_dense(&Atb_F, c);

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
		
#if 0
		double* x_full = (double*)(x->x);
		double* y_full = (double*)(y->x);
		DPRINT(x_full, nvar, " %- .1e");
		DPRINT(y_full, nvar, " %- .1e");
#endif
		cholmod_l_free_dense(&x_F, c);
	}

	cholmod_l_free_dense(&y, c);

	return (x);
}

/*
 * An implementation of the Portugal/Judice/Vicente block-pivoting algorithm for
 * pre-formulated normal equations, with single-row up/down-dates.
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
nnls_normal_block_updown(cholmod_sparse *AtA, cholmod_dense *Atb, int verbose,
   cholmod_common *c)
{
	int nvar = AtA->nrow;
	long F[nvar], G[nvar], H1[nvar], H2[nvar];
	cholmod_dense *x, *y;
	cholmod_factor *L;
	long nF, nG, nH1, nH2, ninf;
	int i, trials, murty_steps, iter;
	clock_t t0, t1;

	/* XXX: make these settable? */
	iter = 3*nvar;		/* Maximum number of iterations */
	trials = MAX_TRIALS;	/* Runs without progress before reverting
				 * to a deterministic algorithm */
	murty_steps = MAX_TRIALS;

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

	/* Drop small entries from the problem, replacing them with eels. */
	//cholmod_l_drop(DBL_EPSILON, AtA, c);

	/* Set up the dual vector */
	for (i = 0; i < nvar; i++)
		((double *)(y->x))[i] = -((double *)(Atb->x))[i];

	AtA->stype = 1;
	t0 = clock();

	/*
	 * Compute a permutation on AtA that approximately minimizes the number
	 * of nonzero elements in its Cholesky decomposition. This permutation
	 * will be used for the remainder of the problem.
	 */
	L = cholmod_l_analyze(AtA, c);

	if (verbose) {
		t1 = clock();
		printf("Analyze[%d]: %.2f s\n",nvar,
		    (double)(t1-t0)/(CLOCKS_PER_SEC)); 
	}

	while (iter-- > 0) {
		/*
		 * Fill H1, H2 with the points that must be moved from the
		 * passive to active and active to passive sets, respectively.
		 */

		nH1 = nH2 = 0;
		for (i = 0; i < nF; i++)
			if (((double *)(x->x))[F[i]] < -KKT_TOL)
				H1[nH1++] = F[i];
		for (i = 0; i < nG; i++)
			if (((double *)(y->x))[G[i]] < -KKT_TOL)
				H2[nH2++] = G[i];

		/*
		 * If all coefficients were positive, then we are done.
		 */

		if (nH1 == 0 && nH2 == 0)
			break;

		if (ninf <= murty_steps)
			trials = -1;

		/*
		 * Check the status of the bulk set switching.
		 * After MAX_TRIALS iterations of Murty's finite method,
		 * revert to block switching.
		 */
		
		if (ninf > murty_steps &&
		    (nH2 + nH1 < ninf || trials < -murty_steps)) {
			if (nH2 + nH1 <= ninf)
			    murty_steps++;
			ninf = nH2 + nH1;
			trials = MAX_TRIALS;
		} else {
			/* Stuck, check if we need to try something else */
			trials--;
			if (verbose)
				printf("\tStuck! trials: %d nH1: %ld nH2: %ld"
				    "\n", trials, nH1, nH2);
			if (trials < 0) {
				/*
				 * Out of luck. Fall back to slow but
				 * guaranteed method (picking the last
				 * infeasible coordinate).
				 */
				
				if (nH2 == 0) {
					goto maxh1;
				} else if (nH1 == 0) {
					goto maxh2;
				} else if (H1[nH1 - 1] > H2[nH2 - 1]) {
					maxh1:
					H1[0] = H1[nH1 - 1];
					nH1 = 1; nH2 = 0;
					if (verbose)
						printf("\tH1: %ld (%e)\n",
						    H1[0],
						    ((double *)(x->x))[H1[0]]);
				} else {
					maxh2:
					H2[0] = H2[nH2 - 1];
					nH2 = 1; nH1 = 0;
					if (verbose)
						printf("\tH2: %ld (%e)\n",
						    H2[0],
						    ((double *)(y->x))[H2[0]]);
				}
			}
		}

		if (verbose)
			printf("Iteration %d Infeasibles: %ld\n",
			    (3*nvar - iter), ninf);

		/* We're about to recompute x, so free it */
		cholmod_l_free_dense(&x, c);

		if (verbose)
		    printf("\tUnconstrained solve for %ld of %d coefficients\n",
		      nF - nH1 + nH2, nvar);

		/*
		 * Next, update the factorization of AtA[:,F][F,:], removing
		 * rows in H1 and adding those in H2. For large updates, the
		 * factorization will be recomputed from scratch; for small
		 * updates, the rows will be added and deleted from the
		 * factorization one-by-one, potentially saving a lot of
		 * computing time.
		 *
		 * F, G, H1, H2, and the associated counts are updated on
		 * exit.
		 */
		L = modify_factor(AtA, L, F, &nF, G, &nG, H1, &nH1, H2, &nH2,
		    verbose, c);
	
		if (verbose) t0 = clock();
		
		/*
		 * Solve the full system, but with the rows of L corresponding
		 * to G set to identity. 
		 */
		x = cholmod_l_solve(CHOLMOD_A, L, Atb, c);

		if (verbose) {
			t1 = clock();
			printf("\tSolve[%d] (%ld free): %.2f s\n",nvar,nF,
			    (double)(t1-t0)/(CLOCKS_PER_SEC));
		}

		/* Explicitly zero x[G] and y[F] */
		for (i = 0; i < nG; i++)
			((double *)(x->x))[G[i]] = 0;
		for (i = 0; i < nF; i++)
			((double *)(y->x))[F[i]] = 0;

		/* Update the constrained part */
		{
			cholmod_sparse *AtA_FG;
			cholmod_dense *Atb_G, *x_F;
		
			/*
			 * NB: We happen to know that AtA is actually
			 * symmetric, but cholmod_submatrix will only
			 * handle unsymmetric ones. The submatrix constructed
			 * this way will also be secretly symmetric.
			 */
			AtA->stype = 0;
			double ones[2] = {1., 0}, mones[2] = {-1., 0};

			if (verbose) t0 = clock();

			Atb_G = cholmod_l_allocate_dense(nG, 1, nG,
			    CHOLMOD_REAL, c);
			for (i = 0; i < nG; i++)
				((double *)(Atb_G->x))[i] =
				    ((double *)(Atb->x))[G[i]];

			x_F = cholmod_l_allocate_dense(nF, 1, nF,
			    CHOLMOD_REAL, c);
			for (i = 0; i < nF; i++)
				((double *)(x_F->x))[i] =
				    ((double *)(x->x))[F[i]];

			AtA_FG = cholmod_l_submatrix(AtA, G, nG, F, nF, 1,1,c);
			cholmod_l_sdmult(AtA_FG, 0, ones, mones, x_F, Atb_G, c);
			
			for (i = 0; i < nG; i++)
				((double *)(y->x))[G[i]] =
				    ((double *)(Atb_G->x))[i];

			cholmod_l_free_dense(&Atb_G, c);
			cholmod_l_free_dense(&x_F, c);
			cholmod_l_free_sparse(&AtA_FG, c);

			if (verbose) {
				t1 = clock();
				printf("\tUpdate y[%ld]: %.2f s\n",
				    nG,(double)(t1-t0)/(CLOCKS_PER_SEC));
			}
			AtA->stype = 1;
			
		}
#if 0
		double* x_full = (double*)(x->x);
		double* y_full = (double*)(y->x);
		DPRINT(x_full, nvar, " %- .1e");
		DPRINT(y_full, nvar, " %- .1e");
#endif
	}

	cholmod_l_free_dense(&y, c);
	cholmod_l_free_factor(&L, c);

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

