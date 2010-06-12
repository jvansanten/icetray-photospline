#include <stdio.h>
//#include <string.h>
#include <math.h>

#include <time.h>
#include <float.h>
#include <stdbool.h>

#include <suitesparse/cholmod.h>
#include "splineutil.h"

static int
intcmp(const void *xa, const void *xb);

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

	cholmod_l_free_factor(&L, c);
	if (n_resolves > 0) cholmod_l_free_dense(&delta_Atb, c);

	return (x);
}

cholmod_sparse* get_column(cholmod_sparse *A, long k, 
    long *iPerm, long *Fset, long nF, cholmod_common *c)
{
	cholmod_sparse *R;
	long *Ap, *Ai, *Anz, *Rp, *Ri;
	double *Ax, *Rx;
	long row;
	int i, j, nz, A_col_nz;

	/* n-by-1, rows not sorted, packed, stype 0 (unsymmetric) */
	R = cholmod_l_allocate_sparse(A->nrow, 1, A->nrow, 
	    false, true, 0, CHOLMOD_REAL, c);

	Ap = (long*)(A->p);   /* pointers to the start of each column */
	Ai = (long*)(A->i);   /* row indices */
	Ax = (double*)(A->x); /* numerical values */
	Anz = (long*)(A->nz); /* length of each column */

	/* get the number of non-zero entries
	 * in column k
	 */
	if (A->packed) A_col_nz = Ap[k+1]-Ap[k];
	else A_col_nz = Anz[k];

	Rp = (long*)(R->p);
	Ri = (long*)(R->i);
	Rx = (double*)(R->x);

	nz = 0;
	/* copy rows from A[F,k] to R */
	for (i = 0; i < nF; i++) {
	/* XXX: O(nF*A_col_nz) brute force search */
		for (j = 0; j < A_col_nz; j++) {
			if (Ai[Ap[k]+j] == Fset[i]) {
				Ri[nz] = Ai[Ap[k]+j];
				Rx[nz] = Ax[Ap[k]+j];
				nz++;	
			}
		} 
	}

	/* symbolically permute the order of the rows */
	if (iPerm) {
		for (i = 0; i < nz; i++) {
			row = Ri[i];
			Ri[i] = iPerm[row];
		}
	}

	/* R is packed */
	Rp[0] = 0;
	Rp[1] = nz;

	return(R);
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

cholmod_factor* modify_factor(cholmod_sparse* A, cholmod_factor *L,
    long *F, long *nF_, long *G, long *nG_, long *H1, long *nH1_,
    long *H2, long *nH2_, bool verbose, cholmod_common *c)
{
	long iPerm[L->n];
	cholmod_sparse *col;
	int i, j, k;
	long nF, nG, nH1, nH2;
	
	nF = *nF_;
	nG = *nG_;
	nH1 = *nH1_;
	nH2 = *nH2_;

	/* Compute the inverse of the fill-reducing permutation */
	for (i = 0; i < L->n; i++)  iPerm[ ((long*)(L->Perm))[i] ] = i;

	/* FIXME: optionally recompute the factorization if it's cheaper */

	/*
	 * Remove elements in H1 from F, and add them to G,
	 * exploiting the fact that H1 elements are in order
	 * relative to their order in F.
	 */
	for (i = 0, j = 0; i < nH1; i++) {
		G[nG++] = H1[i];
		while (F[j] != H1[i]) j++;
		for (k = j+i; k+1 < nF; k++)
			F[k-i] = F[k-i+1];

		/* remove the row from the factorization */
		/* XXX: pass non-zero pattern of row, instead of NULL */
		if (verbose) printf("Deleting row %ld\n", H1[i]);
		cholmod_l_rowdel(iPerm[H1[i]], NULL, L, c);
	}
	nF -= nH1;
	nH1 = 0;

	/* And vice versa */
	for (i = 0, j = 0; i < nH2; i++) {
		F[nF++] = H2[i];
		while (G[j] != H2[i]) j++;
		for (k = j+i; k+1 < nG; k++)
			G[k-i] = G[k-i+1];

		/* extract column H2[i] from A, zeroing any
		 * row not in F and permuting the rows to
		 * match the ordering of the columns in L
		 */
		if (verbose) printf("Adding col %ld\n",H2[i]);
		col = get_column(A, H2[i], iPerm, F, nF, c);
		/* insert the column at the permuted index */
		cholmod_l_rowadd(iPerm[H2[i]], col, L, c);
		cholmod_l_free_sparse(&col, c);
	}
	nG -= nH2;
	nH2 = 0;

	qsort(G, nG, sizeof(G[0]), intcmp);
	qsort(F, nF, sizeof(F[0]), intcmp);

	*nF_  = nF;
	*nG_  = nG;
	*nH1_ = nH1;
	*nH2_ = nH2;
	
	return(L);	
}

