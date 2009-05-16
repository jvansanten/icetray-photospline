#include <stdio.h>

#include "splineutil.h"
#include "glam.h"

/* This isn't a public function */
cholmod_sparse *calc_penalty(int *nsplines, int ndim, int i,
    cholmod_common *c);

int main(void) {
	int nsplines[] = {3,4};
	cholmod_common c;
	cholmod_sparse *P;

	cholmod_l_start(&c);

	P = calc_penalty(nsplines, 2, 1,  &c);
	print_sparse(P, &c);
	cholmod_l_print_sparse(P, "P", &c);
	printf("P stype is %d\n", P->stype);
	
	cholmod_l_finish(&c);
	return 0;
}

