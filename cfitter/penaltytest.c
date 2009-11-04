#include <stdio.h>

#include "splineutil.h"
#include "glam.h"

/* This isn't a public function */
cholmod_sparse *calc_penalty(long *nsplines, double *knots, int ndim, int i,
    int order, int porder, cholmod_common *c);

int main(void) {
	long nsplines[] = {11};
	double knots[] = {1,2,3,4,5,6,7,8,9.5,10,11,11.5,13,14};
	cholmod_common c;
	cholmod_sparse *P;

	cholmod_l_start(&c);

	P = calc_penalty(nsplines, knots, 1, 0, 2, 3, &c);
	print_sparse(P, &c);
	cholmod_l_print_sparse(P, "P", &c);
	printf("P stype is %d\n", P->stype);
	
	cholmod_l_finish(&c);
	return 0;
}

