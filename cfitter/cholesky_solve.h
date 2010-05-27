#ifndef CHOLESKY_SOLVE_H
#define CHOLESKY_SOLVE_H

#include <suitesparse/cholmod.h>

cholmod_dense*
cholesky_solve(cholmod_sparse *AtA, cholmod_dense *Atb, cholmod_common *c, int verbose, int n_resolves);

#endif /* CHOLESKY_SOLVE_H */

