#include <stdio.h>

#include "splineutil.h"
#include "glam.h"

/* This isn't a public function */
cholmod_sparse *calc_penalty(long *nsplines, int ndim, int i, int order,
    cholmod_common *c);

int main(void) {
	long nsplines[] = {6,8};
	cholmod_common c;
	cholmod_sparse *P;

	cholmod_l_start(&c);

	P = calc_penalty(nsplines, 2, 1, 3, &c);
	print_sparse(P, &c);
	cholmod_l_print_sparse(P, "P", &c);
	printf("P stype is %d\n", P->stype);
	
	cholmod_l_finish(&c);
	return 0;
}

