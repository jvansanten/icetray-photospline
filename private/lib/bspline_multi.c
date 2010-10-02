#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__i386__) || defined (__x86_64__)
#include <xmmintrin.h>
#elif defined(__powerpc__)
#include <altivec.h>
#endif

#include "photospline/bspline.h"

#define MAXDIM	    8
#define VECTOR_SIZE 4
typedef float v4sf __attribute__((vector_size(VECTOR_SIZE*sizeof(float))));

#define NVECS MAXDIM/VECTOR_SIZE

#if defined(__i386__) || defined (__x86_64__)
#define v4sf_init(a, b) a = _mm_set1_ps(b)
#elif defined(__powerpc__)
#define v4sf_init(a) a = vec_splats(b)
#else
#define v4sf_init(a, b) { \
	((float *)(&a))[0] = b; \
	((float *)(&a))[1] = b; \
	((float *)(&a))[2] = b; \
	((float *)(&a))[3] = b; \
}
#endif

static int
maxorder(int *order, int ndim)
{
	int i, max = 0;
	
	for (i = 0; i < ndim; i++)
		if (order[i] > max)
			max = order[i];
	
	return (max);
}

static void 
localbasis_multisub(const struct splinetable *table, const int *centers,
    int n, int *restrict pos, unsigned long stride,
    const v4sf **restrict localbasis[table->ndim], v4sf *restrict acc,
    /* On GCC 4.2.1 the stack (so acc_local) is misaligned without this: */
    int gccsucks)
{
	v4sf acc_local[NVECS];
	int i, k;
	
	if (n+1 == table->ndim) {
		/*
		 * If we are at the last recursion level, the weights are
		 * linear in memory, so grab the row-level basis functions
		 * and multiply by the weights. Hopefully the compiler
		 * vector optimizations pick up on this code.
		 */
		
		const int woff = stride + centers[n] - table->order[n];
		const v4sf *row;
		v4sf weight_vec;

		for (k = 0; k <= table->order[n]; k++) {
			row = localbasis[n][k];
			v4sf_init(weight_vec, table->coefficients[woff+k]);
			for (i = 0; i < NVECS; i++)
				acc[i] += row[i]*weight_vec;
		}
	} else {
		for (k = -table->order[n]; k <= 0; k++) {
			/*
			 * If we are not at the last dimension, record where we
			 * are, multiply in the row basis value, and recurse.
			 */

			const v4sf *basis_row = localbasis[n][k+table->order[n]];
			const float *basis_row_ptr =
			    (float*)localbasis[n][k+table->order[n]];
			pos[n] = centers[n] + k;
			for (i = 0; i < NVECS; i++)
				v4sf_init(acc_local[i], 0);
				
			localbasis_multisub(table, centers, n+1, pos,
			    stride + pos[n]*table->strides[n], localbasis,
			    acc_local, 0);
			for (i = 0; i < NVECS; i++) 
				acc[i] += acc_local[i] * basis_row[i];
		}
	}
}

/* Evaluate the spline surface and all its derivatives at x */

void
ndsplineeval_gradient(struct splinetable *table, const double *x,
    const int *centers, double evaluates[table->ndim + 1])
{
	int n, i, j, offset;
	int maxdegree = maxorder(table->order, table->ndim) + 1;
	int nbases = table->ndim + 1;
	int pos[table->ndim];
	v4sf acc[NVECS];
	float valbasis[maxdegree];
	float gradbasis[maxdegree];
	v4sf localbasis[table->ndim][maxdegree][NVECS];
	float *acc_ptr;
	const v4sf *localbasis_rowptr[table->ndim][maxdegree];
	const v4sf **localbasis_ptr[table->ndim];

	if (table->ndim+1 > MAXDIM) {
		fprintf(stderr, "Error: ndsplineeval_gradient() can only "
		    "process up to %d-dimensional tables. Adjust MAXDIM in "
		    "bspline_multi.c to change this.\n", MAXDIM-1);
		exit(1);
	}

		
	for (n = 0; n < table->ndim; n++) {
	
		/* FIXME: compute value and gradient bases in one go */
		bsplvb(table->knots[n], x[n], centers[n],
		    table->order[n] + 1, valbasis);
		bspline_deriv_nonzero(table->knots[n], x[n], centers[n],
		    table->order[n], gradbasis);
		
		for (i = 0; i <= table->order[n]; i++) {
			
			((float*)(localbasis[n][i]))[0] = valbasis[i];
			
			for (j = 1; j < table->ndim+1; j++) {
				if (j == 1+n)
					((float*)(localbasis[n][i]))[j] =
					    gradbasis[i];
				else
					((float*)(localbasis[n][i]))[j] =
					    valbasis[i];
			}
			
			localbasis_rowptr[n][i] = localbasis[n][i];
		}
		
		localbasis_ptr[n] = localbasis_rowptr[n];
	}

	acc_ptr = (float*)acc;

	for (i = 0; i < nbases; i++)
		acc_ptr[i] = 0;

	localbasis_multisub(table, centers, 0, pos, 0, localbasis_ptr, acc, 0);

	for (i = 0; i < nbases; i++)
		evaluates[i] = acc_ptr[i];
}
