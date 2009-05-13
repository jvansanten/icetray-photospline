#include "splineutil.h"

int main(void) {
	cholmod_sparse *a, *b;
	cholmod_common c;

	cholmod_start(&c);

	a = cholmod_read_sparse(stdin, &c);

	printf("Source matrix:\n");
	print_sparse(a, &c);

	b = box(a, a, &c);

	printf("Result:\n");
	print_sparse(b, &c);

	cholmod_finish(&c);

	return (0);
}
