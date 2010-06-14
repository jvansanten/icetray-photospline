#ifndef CHOLESKY_SOLVE_H
#define CHOLESKY_SOLVE_H

#include <suitesparse/cholmod.h>

cholmod_dense*
cholesky_solve(cholmod_sparse *AtA, cholmod_dense *Atb, cholmod_common *c, int verbose, int n_resolves);

cholmod_sparse*
get_column(cholmod_sparse *A, long k, long *iPerm, 
    long *Fset, long nF, cholmod_common *c);

cholmod_factor* 
modify_factor(cholmod_sparse* A, cholmod_factor *L,
    long *F, long *nF, long *G, long *nG, long *H1, long *nH1,
    long *H2, long *nH2, bool verbose, cholmod_common *c);

cholmod_factor* 
recompute_factor(cholmod_sparse *A, cholmod_factor *L, long *iPerm,
    long *F, long nF, cholmod_common *c);

#endif /* CHOLESKY_SOLVE_H */

