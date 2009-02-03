#include <stdio.h>
#include <stdlib.h>

#include "bspline.h"

int main(int argc, char **argv) {
	double knots[] = {3,4,5,6,7,8};
	double x;
	int nknots = sizeof(knots)/sizeof(double);

	x = atof(argv[1]);

	printf("Spline evaluated at %lf:\n",x);
	printf("\t%lf\n",splineeval(knots, nknots, x, 2, -1));

	return 0;
}

	
