#include <stdio.h>
//#include <string.h>
#include <math.h>

#include <time.h>
#include <float.h>
#include <stdbool.h>
#include <assert.h>

#include <suitesparse/cholmod.h>
#include "splineutil.h"
#include "cholesky_solve.h"

static int
intcmp(const void *xa, const void *xb);

cholmod_dense *
cholesky_solve(cholmod_sparse *AtA, cholmod_dense *Atb, cholmod_common *c,
     int verbose, int n_resolves)
{
	int i,j,nvar;
	double ones[2] = {1., 0}, mones[2] = {-1., 0};
	double sum_delta_b;
	cholmod_dense *x, *delta_x, *delta_Atb;
	cholmod_factor *L;
	clock_t t0, t1;

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
	R = cholmod_l_allocate_sparse(A->nrow, 1, nF, 
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

cholmod_factor* 
modify_factor(cholmod_sparse *A, cholmod_factor *L,
    long *F, long *nF_, long *G, long *nG_, long *H1, long *nH1_,
    long *H2, long *nH2_, bool verbose, cholmod_common *c)
{
	double changed;
	bool update;
	long nF, nG, nH1, nH2;

	nF = *nF_;
	nG = *nG_;
	nH1 = *nH1_;
	nH2 = *nH2_;

	update = true;
	changed = 1.0;
	/*
	 * If available, use the flop counts from the most recent call to
	 * cholmod_analyze or cholmod_rowadd/rowdel to estimate the expense
	 * of the job.
	 *
	 * XXX: assumptions here are:
	 * 1. we can get ~ 8x speedup from 16 threads in GotoBLAS. 
	 * 2. the size of F is roughly constant: the next factorization will
	 *    be about as expensive as the previous one.
	 */
	if ((c->fl > 0) && (c->modfl > 0)) {
		update = ((c->fl)/8.0 > ((double)(nH1 + nH2))*(c->modfl));
		if (verbose)
			printf("\tFactor flops: %.0lf Mod flops: %.0lf\n",
			    c->fl, c->modfl);
	} else {
		/* 
	 	 * This is th Portugal/Judice/Vincente heuristic for updating
		 * a QR factorization, and not necessarily related to the work
		 * required by CHOLMOD.
		 */
		changed = ((double)(nH1 + nH2))/((double)(nF - nH1 + nH2));
		if ((nF == 0) || (changed > 0.2)) update = false;
	}

	/* short-circuit for rank-1 updates */
	if (nH1 + nH2 == 1) update = true;

	if ((!update) && verbose)
		printf("\tRecomputing factorization from scratch "
		    "(F[%ld], G[%ld], H1[%ld], H2[%ld], changed=%.2lf)\n",
		    nF, nG, nH1, nH2, changed);

	return(modify_factor_p(A, L, F, nF_, G, nG_, H1, nH1_, H2, nH2_,
	    update, verbose, c));

}

cholmod_factor* 
modify_factor_p(cholmod_sparse *A, cholmod_factor *L,
    long *F, long *nF_, long *G, long *nG_, long *H1, long *nH1_,
    long *H2, long *nH2_, bool update, bool verbose, cholmod_common *c)
{
	long iPerm[L->n];
	cholmod_sparse *col;
	int i, j, k;
	long nF, nG, nH1, nH2;
	double changed;
	clock_t t0, t1;
	
	nF = *nF_;
	nG = *nG_;
	nH1 = *nH1_;
	nH2 = *nH2_;

	/* Compute the inverse of the fill-reducing permutation */
	for (i = 0; i < L->n; i++)  iPerm[ ((long*)(L->Perm))[i] ] = i;


	t0 = clock();

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

		if (update) {
			/* remove the row from the factorization */
			/* XXX: pass non-zero pattern of row, instead of NULL */
			cholmod_l_rowdel(iPerm[H1[i]], NULL, L, c);
		}
	}
	if (verbose && update && (nH1 > 0)) {
		t1 = clock();
		printf("\tDelete %ld rows: %.2f s\n",nH1,
		    (double)(t1-t0)/(CLOCKS_PER_SEC));
	}
	nF -= nH1;
	nH1 = 0;

	t0 = clock();

	/* And vice versa */
	for (i = 0, j = 0; i < nH2; i++) {
		F[nF++] = H2[i];
		while (G[j] != H2[i]) j++;
		for (k = j+i; k+1 < nG; k++)
			G[k-i] = G[k-i+1];

		if (update) {
			/* 
			 * Extract column H2[i] from A, zeroing any
			 * row not in F and permuting the rows to
			 * match the ordering of the columns in L.
			 */
			col = get_column(A, H2[i], iPerm, F, nF, c);
			/* insert the column at the permuted index */
			cholmod_l_rowadd(iPerm[H2[i]], col, L, c);
			cholmod_l_free_sparse(&col, c);
		}
	}
	if (verbose && update && (nH2 > 0)) {
		t1 = clock();
		printf("\tAdd %ld rows: %.2f s\n",nH2,
		    (double)(t1-t0)/(CLOCKS_PER_SEC));
	}
	nG -= nH2;
	nH2 = 0;

	qsort(G, nG, sizeof(G[0]), intcmp);
	qsort(F, nF, sizeof(F[0]), intcmp);

	/*
	 * If no update was performed, recompute the
	 * factorization of A[:,F][F,:] from scratch.
	 */
	if (!update) {
		t0 = clock();
		L = recompute_factor(A, L, iPerm, F, nF, c);
		if (verbose) {
			t1 = clock();
			printf("\tFactorize[%ld]: %.2f s\n",nF,
			    (double)(t1-t0)/(CLOCKS_PER_SEC));
		}
	}

	*nF_  = nF;
	*nG_  = nG;
	*nH1_ = nH1;
	*nH2_ = nH2;
	
	return(L);	
}

cholmod_factor* 
recompute_factor(cholmod_sparse *A, cholmod_factor *L, long *iPerm,
    long *F, long nF, cholmod_common *c)
{
	long iF[L->n];
	long FPerm[nF];
	long Lrows[nF]; /* mapping from permuted Lprime to permuted L */
	long *LPerm, *LColCount, *L_FColCount;
	long *Li, *Lp, *Lnz, *Lnext, *L_Fi, *L_Fp, *L_Fnz;
	double *Lx, *L_Fx;
	cholmod_sparse *A_F;
	cholmod_factor *L_F;
	int i, j, nz, nFPerm, common_nmethods, Astype;
	bool common_postorder;

	/* 
	 * Scortched earth: free the unneeded numeric 
	 * bits of L, should they exist.
	 */
	cholmod_l_change_factor(CHOLMOD_PATTERN, 
	    false, /* to LL' */
	    false, /* to supernodal */
	    false, /* to packed */
	    false, /* to monotonic */
	    L, c);

	LPerm = (long*)(L->Perm);
	LColCount = (long*)(L->ColCount);

	/* build an inverse mapping for F */
	for (i = 0; i < L->n; i++) iF[i] = -1;
	for (i = 0; i < nF; i++) iF[F[i]] = i;

	/* Permute F to match L->Perm */
	nFPerm = 0;
	if (iPerm) {
		/* calculate the permuation as it applies to subset F */
		for (i = 0; i < L->n; i++) {
			if (iF[LPerm[i]] >= 0) {
				FPerm[nFPerm] = iF[LPerm[i]];
				/* 
				 * This row of L_F (corrensponding to F[j]) 
				 * corresponds to the position of F[j] in the
				 * fill-reducing permuation.
				 */
				Lrows[nFPerm] = iPerm[F[FPerm[nFPerm]]];
				nFPerm++;
			}
		}
		
	} else {
		/* mostly a no-op in natural ordering */
		for (i = 0; i < nF; i++) {
			FPerm[nFPerm] = i;
			Lrows[nFPerm] = F[i];
			nFPerm++;
		}
	}

	Astype = A->stype;
	assert( Astype != 0);
	A->stype = 0; /* cholmod_submatrix doesn't like symmetric matrices */
	A_F = cholmod_l_submatrix(A, F, nF, F, nF, 1, 1, c);
	A->stype = Astype;
	A_F->stype = Astype;
	A_F->stype = 1;

	/* 
	 * Since we intend to copy the results of the submatrix factorization
	 * in to L, we need to force cholmod_analyze to use exactly the same
	 * permutation previously calculated for L, without postordering.
	 */
	common_nmethods = c->nmethods;
	common_postorder = c->postorder;
	c->nmethods = 1;
	c->postorder = false;

	/* calculate pattern of L in the given permuation */
	L_F = cholmod_l_analyze_p(A_F, FPerm, NULL, -1, c);
	L_FColCount = (long*)(L_F->ColCount);

	/* restore alternate orderings */
	c->nmethods = common_nmethods;
	c->postorder = common_postorder;

#ifndef NDEBUG
	/* sanity check: did we get the permutation we asked for? */
	long* L_FPerm = (long*)(L_F->Perm);
	for (i = 0; i < nF; i++) assert( FPerm[i] == L_FPerm[i] );
#endif

	cholmod_l_factorize(A_F, L_F, c);

	cholmod_l_free_sparse(&A_F, c);

	/* We really, really, really don't want L in supernodal form. */
	cholmod_l_change_factor(CHOLMOD_REAL, 
	    false, /* to LL' */
	    false, /* to supernodal */
	    false, /* to packed */
	    false, /* to monotonic */
	    L_F, c);

	/* 
	 * XXX: swizzle the nonzero pattern of L to match
	 * L'. Does this break anything down the line?
	 */
	for (i = 0; i < L->n; i++) LColCount[i] = 1;
	for (i = 0; i < nF; i++) LColCount[Lrows[i]] = L_FColCount[i];

	/* 
	 * Restore numeric bits of L, initialized to identity 
	 * and with the new nonzero pattern.
	 */
	cholmod_l_change_factor(CHOLMOD_REAL, 
	    false, /* to LL' */
	    false, /* to supernodal */
	    false, /* to packed */
	    false, /* to monotonic */
	    L, c);

	/* 
	 * Now that we have simplicial numeric factorizations
	 * of both L_F and L (currently identity), we can get 
	 * the appropriate pointers to the numerical bits of both.
	 */
	Lp =    (long*)(L->p);
	Lnz =   (long*)(L->nz);
	Li =    (long*)(L->i);
	Lnext = (long*)(L->next);
	Lx =    (double*)(L->x);
	L_Fp =  (long*)(L_F->p);
	L_Fnz = (long*)(L_F->nz);
	L_Fi =  (long*)(L_F->i);
	L_Fx =  (double*)(L_F->x);

	/* 
	 * Copy the numeric factorization from L_F to 
	 * the corresponding columns and rows of L.
	 */
	for (i = 0; i < nF; i++) {
		/* 
		 * NB: cholmod_factor is never explicitly packed
		 * (unlike cholmod_sparse, L->nz is always meaningful).
		 */
		if (L_Fnz) nz = L_Fnz[i];
		else nz = L_Fp[i+1] - L_Fp[i];

		if ( Lp[Lnext[Lrows[i]]] - Lp[Lrows[i]] < nz ) {
			cholmod_l_reallocate_column(Lrows[i], nz, L, c);
			printf("L->nz[%ld] <= %ld, L_F->nz[%d] = %ld\n", 
		    	    Lrows[i], Lp[Lnext[Lrows[i]]] - Lp[Lrows[i]],
			    i, L_Fnz[i]);
		}

		assert( Lp[Lnext[Lrows[i]]] - Lp[Lrows[i]] >= nz );

		for (j = 0; j < nz; j++) {
			/*
			 * Adjust index of each row to match
			 * its position in L.
			 */
			Li[Lp[Lrows[i]]+j] = Lrows[L_Fi[L_Fp[i]+j]] ;	
			/*
			 * Numerical values can just be copied
			 * to the appropriate column.
			 */
			Lx[Lp[Lrows[i]]+j] = L_Fx[L_Fp[i]+j];	
		}
		Lnz[Lrows[i]] = nz;	
	}

	/* we're done with L_F */
	cholmod_l_free_factor(&L_F, c);

	return(L);
}

