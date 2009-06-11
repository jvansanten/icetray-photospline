#include <stdio.h>
#include <string.h>
#include <math.h>

#include <suitesparse/SuiteSparseQR_C.h>

#include "splineutil.h"
#include "splinetable.h"
#include "glam.h"

static cholmod_sparse *flatten_ndarray_to_sparse(struct ndsparse *array,
    size_t nrow, size_t ncol, cholmod_common *c);
cholmod_sparse *calc_penalty(long *nsplines, int ndim, int i,
    cholmod_common *c);
void print_ndsparse_py(struct ndsparse *a);

#define max(a,b) ((a > b) ? a : b)

void
glamfit(struct ndsparse *data, double *weights, double **coords,
    struct splinetable *out, double smooth, int order, int verbose,
    cholmod_common *c)
{
	cholmod_sparse **bases, **boxedbases;
	cholmod_sparse *penalty, *penalty_chunk;
	cholmod_dense *coefficients, *Rdens;
	struct ndsparse F, R;
	cholmod_sparse *Fmat, *Rmat, *fitmat;
	double scale1[2] = {1.0, 0.0}, scale2[2] = {1.0, 0.0};
	size_t sidelen;
	long *nsplines;
	long i, j, n;

	nsplines = calloc(data->ndim,sizeof(long));

	/*
	 * We will start by precomputing some things before we don't need to
	 * worry about memory pressure on temporaries.
	 */

	/* Figure out the total number of spline nodes */
	sidelen = 1;
	for (i = 0; i < out->ndim; i++) {
		nsplines[i] = out->nknots[i] - order - 1;
		sidelen *= nsplines[i];
	}

	/*
	 * Compute the spline basis matrices, and pre-box them.
	 */

	if (verbose)
		printf("Calculating spline basis...\n");

	bases = calloc(out->ndim, sizeof(cholmod_sparse *));
	boxedbases = calloc(out->ndim, sizeof(cholmod_sparse *));

	for (i = 0; i < out->ndim; i++) {
		bases[i] = bsplinebasis(out->knots[i], out->nknots[i], coords[i], 
		    data->ranges[i], order, c);
		boxedbases[i] = box(bases[i], bases[i], c);
	}

	/* 
	 * Now we want to get the penalty matrix.
	 */

	if (verbose)
		printf("Calculating penalty matrix...\n");

	penalty = cholmod_l_spzeros(sidelen, sidelen, 1, CHOLMOD_REAL, c);
	for (i = 0; i < out->ndim; i++) {
		cholmod_sparse *penalty_tmp;

		penalty_chunk = calc_penalty(nsplines, out->ndim, i, c);
		penalty_tmp = penalty;

		/* Add each chunk to the big matrix, scaling by smooth */
		scale2[0] = smooth;
		penalty = cholmod_l_add(penalty, penalty_chunk, scale1, scale2,
		    1, 0, c);

		cholmod_l_free_sparse(&penalty_tmp, c);
		cholmod_l_free_sparse(&penalty_chunk, c);
	}

	
	/*
	 * Begin iterating.
	 *
	 * XXX: clamped to 1 iteration
	 */

	if (verbose)
		printf("Reticulating splines...\n");

	for (n = 0; n < 1; n++) {	/* XXX: actual fit iteration unimplemented */
		/*
		 * Initialize F and R. 
		 * F = weights
		 * R = weights * data
		 */
		R.rows = F.rows = data->rows;
		R.ndim = F.ndim = data->ndim;
		F.x = malloc(data->rows * sizeof(double));
		F.i = malloc(2*data->ndim * sizeof(int *));
		F.ranges = malloc(2*data->ndim * sizeof(int));
		R.x = malloc(data->rows * sizeof(double));
		R.i = malloc(data->ndim * sizeof(int *));
		R.ranges = malloc(data->ndim * sizeof(int));

		for (i = 0; i < data->ndim; i++) {
			F.i[i] = malloc(data->rows * sizeof(int));
			R.i[i] = malloc(data->rows * sizeof(int));
		}
	
		memcpy(R.x, weights, data->rows * sizeof(double));
		memcpy(F.x, weights, data->rows * sizeof(double));

		for (i = 0; i < data->ndim; i++) {
			memcpy(R.i[i], data->i[i], data->rows * sizeof(int));
			memcpy(F.i[i], data->i[i], data->rows * sizeof(int));
			R.ranges[i] = F.ranges[i] = data->ranges[i];
		}

		for (i = 0; i < data->rows; i++)
			R.x[i] *= data->x[i];

		/*
		 * Convolve F and R with the basis matrices
		 */

		if (verbose)
			printf("\tConvolving bases...\n");

		for (i = 0; i < data->ndim; i++) {
			if (verbose)
				printf("\t\tConvolving dimension %ld\n",i);

			slicemultiply(&F, boxedbases[i], i, c);
			slicemultiply(&R, bases[i], i, c);
		}

		/* Now flatten R into a matrix */

		if (verbose)
			printf("\tFlattening residuals matrix...\n");

		Rmat = flatten_ndarray_to_sparse(&R, sidelen, 1, c);

		for (i = 0; i < R.ndim; i++)
			free(R.i[i]);
		free(R.x); free(R.i); free(R.ranges);

		/* XXX: We reshape, transpose, and then flatten F, which is
		 * potentially memory hungry. This can probably be done in one
		 * step. */

		/*
		 * F now has the number of splines squared as each axis dimension.
		 * We now want to double the dimensionality of F so that it is
		 * n1xn1xn2xn2x... instead of (n1xn1)x(n2xn2)x...
		 */

		if (verbose)
			printf("Transforming fit array...\n");

		/* Fill the ranges array in-place by starting at the back, and make
		 * use of 3/2 = 2/2 = 1 to get the pairs. While here, rearrange the
		 * index columns using the same logic. */
		F.ndim *= 2;
		for (i = F.ndim-1; i >= 0; i--) {
			F.ranges[i] = sqrt(F.ranges[i/2]);
			if (i % 2 == 0)
				F.i[i] = F.i[i/2];
			else
				F.i[i] = malloc(F.rows * sizeof(int));
		}

		/* Now figure out each point's new coordinates */
		for (i = 0; i < F.rows; i++) {
			for (j = 0; j < F.ndim; j += 2) {
				F.i[j+1][i] = F.i[j][i] % F.ranges[j];
				F.i[j][i] = F.i[j][i] / F.ranges[j];
			}
		}

		/* Transpose F so that the even-numbered axes come first */
		{
			int **oldi, *oldranges;

			oldi = F.i;
			oldranges = F.ranges;
			F.i = malloc(F.ndim * sizeof(int *));
			F.ranges = malloc(F.ndim * sizeof(int));
			for (i = 0; i < F.ndim; i++) {
				if (i % 2 == 0) {
					F.i[i/2] = oldi[i];
					F.ranges[i/2] = oldranges[i];
				} else {
					F.i[F.ndim/2 + i/2] = oldi[i];
					F.ranges[F.ndim/2 + i/2] = oldranges[i];
				}
			}
			free(oldi);
			free(oldranges);
		}

		/* Now flatten F */

		Fmat = flatten_ndarray_to_sparse(&F, sidelen, sidelen, c);
		for (i = 0; i < F.ndim; i++)
			free(F.i[i]);
		free(F.x); free(F.i); free(F.ranges);
	
		/* XXX: optimization possibilities ended */

		scale2[0] = 1.0;
		fitmat = cholmod_l_add(Fmat, penalty, scale1, scale2, 1, 0, c);
		cholmod_l_free_sparse(&Fmat, c); /* we don't need Fmat anymore */

		Rdens = cholmod_l_sparse_to_dense(Rmat, c);
		cholmod_l_free_sparse(&Rmat, c); /* we don't need Fmat anymore */

		/*
		 * Now, we can solve the linear system 
		 */

		if (verbose)
			printf("Computing iteration %ld least square solution...\n",
			    n+1);
	
		coefficients = SuiteSparseQR_C_backslash_default(fitmat, Rdens, c);
		cholmod_l_free_sparse(&fitmat, c);
		if (coefficients == NULL)
			printf("Solution FAILED\n");
	}
	
	/* Clean up detritus */

	if (verbose)
		printf("Done: cleaning up\n");

	cholmod_l_free_sparse(&Fmat, c);
	cholmod_l_free_dense(&Rdens, c);
	cholmod_l_free_sparse(&penalty, c);

	for (i = 0; i < out->ndim; i++) {
		cholmod_l_free_sparse(&bases[i], c);
		cholmod_l_free_sparse(&boxedbases[i], c);
	}
	free(bases); free(boxedbases);
	free(nsplines);

	/* Copy out the coefficients */

	out->coefficients = malloc(coefficients->nrow * coefficients->ncol *
	    sizeof(double));
	memcpy(out->coefficients, coefficients->x,
	    coefficients->nrow * coefficients->ncol * sizeof(double));
	out->naxes = malloc(out->ndim * sizeof(long));
	for (i = 0; i < out->ndim; i++)
		out->naxes[i] = out->nknots[i] - order - 1;

	/* Free our last matrix */
	cholmod_l_free_dense(&coefficients, c);
}

static cholmod_sparse *
flatten_ndarray_to_sparse(struct ndsparse *array, size_t nrow, size_t ncol,
    cholmod_common *c)
{
	cholmod_triplet *trip;
	cholmod_sparse *sparse;
	long moduli[array->ndim];
	long i, j, k;

	trip = cholmod_l_allocate_triplet(nrow, ncol, array->rows, 0,
	    CHOLMOD_REAL, c);

	moduli[array->ndim-1] = 1;
	for (i = array->ndim-2; i >= 0; i--)
		moduli[i] = moduli[i+1]*array->ranges[i+1];

	for (i = 0; i < array->rows; i++) {
		k = 0;
		for (j = 0; j < array->ndim; j++)
			k += array->i[j][i]*moduli[j];

		((long *)(trip->j))[i] = k % ncol;
		((long *)(trip->i))[i] = k / ncol;
		((double *)(trip->x))[i] = array->x[i];
	}
	trip->nnz = array->rows;

	sparse = cholmod_l_triplet_to_sparse(trip, trip->nnz, c);
	cholmod_l_free_triplet(&trip, c);

	return (sparse);
}

static cholmod_sparse *
offdiagonal(int rows, int cols, int diag, cholmod_common *c)
{
	int row, col;
	cholmod_triplet *trip;
	cholmod_sparse *sparse;
	
	if (diag <= 0) {
		row = diag;
		col = 0;
	} else {
		col = diag;
		row = 0;
	}

	trip = cholmod_l_allocate_triplet(rows, cols,
	    max(rows,cols) - abs(diag), 0, CHOLMOD_REAL, c);

	while (row < rows && col < cols) {
		((long *)(trip->i))[trip->nnz] = row;
		((long *)(trip->j))[trip->nnz] = col;
		((double *)(trip->x))[trip->nnz] = 1.0;

		trip->nnz++; row++; col++;
	}

	sparse = cholmod_l_triplet_to_sparse(trip, trip->nnz, c);
	cholmod_l_free_triplet(&trip, c);

	return sparse;
}

cholmod_sparse *
calc_penalty(long *nsplines, int ndim, int dim, cholmod_common *c)
{
	cholmod_sparse *finitediff, *fd_trans, *DtD;
	cholmod_sparse *tmp, *tmp2, *result;
	double scale1[2] = {1.0, 0.0}, scale2[2] = {1.0, 0.0};
	long i;

	/* First, we will compute the finite difference matrix,
	 * which looks like this:
	 * 
	 * 1 -2  1  0 0 ...
	 * 0  1 -2  1 0 ...
	 * 0  0  1 -2 1 ...
	 */

	tmp = cholmod_l_speye(nsplines[dim] - 2, nsplines[dim],
	    CHOLMOD_REAL, c);
	tmp2 = offdiagonal(nsplines[dim] - 2, nsplines[dim],
	    1, c);
	scale2[0] = -2.0;
	finitediff = cholmod_l_add(tmp, tmp2, scale1, scale2, 1,
	    0, c);
	cholmod_l_free_sparse(&tmp2, c);
	cholmod_l_free_sparse(&tmp, c);

	scale2[0] = 1.0;
	tmp2 = offdiagonal(nsplines[dim] - 2, nsplines[dim],
	    2, c);
	tmp = cholmod_l_add(finitediff, tmp2, scale1, scale2, 1,
	    0, c);
	cholmod_l_free_sparse(&tmp2, c);
	cholmod_l_free_sparse(&finitediff, c);
	finitediff = tmp;

	/*
	 * Now we want DtD, which is the transpose of finitediff
	 * multiplied by finitediff
	 */

	fd_trans = cholmod_l_transpose(finitediff, 1, c);
	DtD = cholmod_l_ssmult(fd_trans, finitediff,
	    1 /* DtD is symmetric */, 1, 0, c);
	cholmod_l_free_sparse(&finitediff, c);
	cholmod_l_free_sparse(&fd_trans, c);

	/* Next take kronecker products to form the full P */

	tmp = NULL;
	result = NULL;
	for (i = 0; i < ndim; i++) {
		tmp2 = (i == dim) ? DtD : cholmod_l_speye(
		    nsplines[i], nsplines[i],
		    CHOLMOD_REAL, c);
		tmp2->stype = 1; /* The identity matrix is always symmetric. */

		if (result == NULL) {
			result = tmp2;
			continue;
		}

		tmp = kronecker_product(result, tmp2, c);
		cholmod_l_free_sparse(&result, c);
		cholmod_l_free_sparse(&tmp2, c);

		result = tmp;
	}

	return result;
}

