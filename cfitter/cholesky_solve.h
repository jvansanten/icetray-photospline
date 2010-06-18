#ifndef CHOLESKY_SOLVE_H
#define CHOLESKY_SOLVE_H

#include <suitesparse/cholmod.h>
#include <stdbool.h>

cholmod_dense *
cholesky_solve(cholmod_sparse *AtA, cholmod_dense *Atb, cholmod_common *c, int verbose, int n_resolves);

cholmod_sparse *
get_column(cholmod_sparse *A, long k, long *iPerm, 
    long *Fset, long nF, cholmod_common *c);

cholmod_sparse *
submatrix_symm(cholmod_sparse *A, long *rows, long nrows,
    long *cols, long ncols, cholmod_common *c);

cholmod_factor * 
modify_factor(cholmod_sparse* A, cholmod_factor *L,
    long *F, long *nF, long *G, long *nG, long *H1, long *nH1,
    long *H2, long *nH2, int verbose, cholmod_common *c);

cholmod_factor * 
modify_factor_p(cholmod_sparse *A, cholmod_factor *L,
    long *F, long *nF_, long *G, long *nG_, long *H1, long *nH1_,
    long *H2, long *nH2_, bool update, bool verbose, cholmod_common *c);

cholmod_factor * 
recompute_factor(cholmod_sparse *A, cholmod_factor *L, long *iPerm,
    long *F, long nF, cholmod_common *c);

#define FACTOR_INFO(L) printf(#L " is_ll: %d, is_super: %d, is_monotonic: %d, xtype: %d, ordering: %d\n",L->is_ll,L->is_super,L->is_monotonic,L->xtype,L->ordering)

#define DPRINT(vec,size,format)\
	{\
	int i;\
	printf(#vec":");\
	for (i = 0; i < size; i++) printf(format,vec[i]);\
	printf("\n");\
	}

#endif /* CHOLESKY_SOLVE_H */

