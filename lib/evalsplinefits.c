#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "splinetable.h"
#include "bspline.h"

void usage() {
	fprintf(stderr,"evalsplinefits <path> x1 x2 x3 ...\n");
	exit(1);
}

int main(int argc, char **argv) {
	struct splinetable table;
	int i;
	double *x;
	int *centers;

	if (argc < 2)
		usage();

	readsplinefitstable(argv[1], &table);

	if (argc < 2+table.ndim)
		usage();

	x = malloc(sizeof(double)*table.ndim);
	centers = malloc(sizeof(double)*table.ndim);
	for (i = 0; i < table.ndim; i++) 
		sscanf(argv[i+2],"%lf",&x[i]);
	tablesearchcenters(&table, x, centers);

	printf("NDim: %d\n",table.ndim);
	printf("Order: %d\n",table.order);

	printf("Value: %lf\n",ndsplineeval(&table, x, centers));

	return 0;
}
