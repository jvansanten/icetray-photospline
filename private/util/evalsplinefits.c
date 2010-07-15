#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "splinetable.h"
#include "bspline.h"

static void usage() {
	fprintf(stderr,"evalsplinefits <path> x1 x2 x3 ...\n");
	exit(1);
}

#define TIMING

int main(int argc, char **argv) {
	struct splinetable table;
	int i;
	double *x, value;
	int *centers;
	struct timeval tp1, tp2;

	if (argc < 2)
		usage();

	gettimeofday(&tp1, NULL);
	readsplinefitstable(argv[1], &table);
	gettimeofday(&tp2, NULL);

    #ifdef TIMING
	if (tp2.tv_usec < tp1.tv_usec)
		tp2.tv_usec += 1e6;
	printf("Time to open table: %ld.%06ld seconds\n", 
	    tp2.tv_sec - tp1.tv_sec, tp2.tv_usec - tp1.tv_usec);
    #endif

	if (argc < 2+table.ndim)
		usage();

	x = malloc(sizeof(double)*table.ndim);
	centers = malloc(sizeof(double)*table.ndim);
	for (i = 0; i < table.ndim; i++) 
		sscanf(argv[i+2],"%lf",&x[i]);
	gettimeofday(&tp1, NULL);
	tablesearchcenters(&table, x, centers);
	value = ndsplineeval(&table, x, centers, 0);
	gettimeofday(&tp2, NULL);

    #ifdef TIMING
	if (tp2.tv_usec < tp1.tv_usec)
		tp2.tv_usec += 1e6;
	printf("Cold cache evaluation time: %ld microseconds\n",
	    tp2.tv_usec - tp1.tv_usec);

	gettimeofday(&tp1, NULL);
	value = ndsplineeval(&table, x, centers, 0);
	gettimeofday(&tp2, NULL);

	if (tp2.tv_usec < tp1.tv_usec)
		tp2.tv_usec += 1e6;
	printf("Warm cache evaluation time: %ld microseconds\n",
	    tp2.tv_usec - tp1.tv_usec);
    #endif

	printf("NDim: %d\n",table.ndim);
	printf("Order: %d\n",table.order[0]);

	printf("Value: %lf\n",value);
	printf("e^(value): %e\n",exp(value));

	return 0;
}
