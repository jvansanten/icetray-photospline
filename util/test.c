#include <stdio.h>
#include <stdlib.h>

#include "bspline.h"

int main(int argc, char **argv) {
	double knots[] = {2,3,4,5,6,7,8,9,10};
	double weights[] = {2,3,4,5,6,7};
	double x;
	int nknots = sizeof(knots)/sizeof(double);
	int centers[] = {3};
	double *gra = knots;
	double **gra2 = &gra;

	x = atof(argv[1]);

	printf("Spline evaluated at %lf:\n",x);
	printf("\t%lf\n",splineeval(knots, weights, nknots, x, 2, centers[0]));
	printf("\t%lf\n",ndsplineeval(gra2, weights, 1, &nknots, &x, 2, centers));

	return 0;
}

	
