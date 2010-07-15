#include "splineutil.h"

int main(void) {
	cholmod_sparse *a, *b;
	cholmod_common c;

	cholmod_l_start(&c);

	a = cholmod_l_read_sparse(stdin, &c);

	printf("Source matrix:\n");
	print_sparse(a, &c);

	b = kronecker_product(a, a, &c);

	printf("Result:\n");
	print_sparse(b, &c);

	cholmod_l_finish(&c);

	return (0);
}
