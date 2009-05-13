#include "splineutil.h"

/* Compare result to:
 * import glam
 * import numpy
 * 
 * b = numpy.reshape(range(0,27),(3,3,3))
 * a = 2*numpy.eye(3)
 * a[0,0] = 1
 * a[0,2] = 1
 * glam.rho(a,b,2)
 */

void printndsparse(struct ndsparse *a) {
	double x[27];
	int i;

	for (i = 0; i < 27; i++)
		x[i] = 0;
	for (i = 0; i < 27; i++)
		x[a->i[2][i] + a->i[1][i]*3 + a->i[0][i]*9] = a->x[i];
	for (i = 0; i < 27; i++) {
		if (i % 3 == 0 && i != 0) printf("\n");
		if (i % 9 == 0 && i != 0) printf("\n");
		printf("%lf\t",x[i]);
	}
	printf("\n");
}


int main(void) {
	struct ndsparse b;
	cholmod_sparse *a;
	cholmod_common c;
	int i;

	cholmod_start(&c);

	b.rows = 27;
	b.ndim = 3;
	b.i = malloc(sizeof(int *)*3);
	b.x = malloc(sizeof(double)*27);
	b.ranges = malloc(sizeof(int)*3);
	for (i = 0; i < 3; i++) {
		b.i[i] = malloc(sizeof(int)*27);
		b.ranges[i] = 3;
	}
	for (i = 0; i < 27; i++) {
		b.x[i] = i;
		b.i[0][i] = i / 9;
		b.i[1][i] = (i / 3) % 3;
		b.i[2][i] = i % 3;
	}

	printf("Input array:\n");
	printndsparse(&b);
	
	a = cholmod_read_sparse(stdin, &c);

	printf("Multiplying on dimension %d by:\n",2);
	print_sparse(a, &c);

	slicemultiply(&b, a, 2, &c);

	printf("Result:\n");
	printndsparse(&b);

	cholmod_finish(&c);

	return (0);
}
