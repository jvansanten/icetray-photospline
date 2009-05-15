#include <stdio.h>
#include <string.h>
#include <math.h>

#include <suitesparse/SuiteSparseQR_C.h>

#include "splineutil.h"
#include "splinetable.h"
#include "glam.h"

static cholmod_sparse *flatten_ndarray_to_sparse(struct ndsparse *array,
    size_t nrow, size_t ncol, cholmod_common *c);
cholmod_sparse *calc_penalty(int *nsplines, int ndim, int i,
    cholmod_common *c);

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
	int *moduli, *nsplines;
	int i, j, n;

	moduli = calloc(2*data->ndim,sizeof(int));
	nsplines = calloc(data->ndim,sizeof(int));

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

	penalty = cholmod_spzeros(sidelen, sidelen, 1, CHOLMOD_REAL, c);
	for (i = 0; i < out->ndim; i++) {
		cholmod_sparse *penalty_tmp;

		penalty_chunk = calc_penalty(nsplines, out->ndim, i, c);
		penalty_tmp = penalty;

		/* Add each chunk to the big matrix, scaling by smooth */
		scale2[0] = smooth;
		penalty = cholmod_add(penalty, penalty_chunk, scale1, scale2, 1, 0, c);

		cholmod_free_sparse(&penalty_tmp, c);
		cholmod_free_sparse(&penalty_chunk, c);
	}

	
	/*
	 * Begin iterating.
	 *
	 * XXX: clamped to 1 iteration
	 */

	if (verbose)
		printf("Reticulating splines...\n");

	/*
	 * Initialize F and R
	 */

	R.rows = F.rows = data->rows;
	R.ndim = F.ndim = data->ndim;
	F.x = malloc(data->rows * sizeof(double));
	F.i = malloc(2*data->ndim * sizeof(int *));
	F.ranges = malloc(data->ndim * sizeof(int));
	R.x = malloc(data->rows * sizeof(double));
	R.i = malloc(data->ndim * sizeof(int *));
	R.ranges = malloc(data->ndim * sizeof(int));

	for (i = 0; i < data->ndim; i++) {
		F.i[i] = malloc(data->rows * sizeof(int));
		R.i[i] = malloc(data->rows * sizeof(int));
	}

	for (n = 0; n < 1; n++) {	/* XXX: actual fit iteration unimplemented */
		/*
		 * Initialize F and R. 
		 * F = weights
		 * R = weights * data
		 */
	
		memcpy(R.x, weights, data->rows * sizeof(double));
		memcpy(F.x, weights, data->rows * sizeof(double));

		for (i = 0; i < data->ndim; i++) {
			memcpy(R.i[i], data->i[i], data->rows * sizeof(int));
			memcpy(F.i[i], data->i[i], data->rows * sizeof(int));
			R.ranges[i] = F.ranges[i] = data->ranges[i];
		}

		for (i = 0; i < data->rows; i++)
			R.x[i] *= weights[i];

		/*
		 * Convolve F and R with the basis matrices
		 */

		if (verbose)
			printf("\tConvolving bases...\n");

		for (i = 0; i < data->ndim; i++) {
			if (verbose)
				printf("\t\tConvolving dimension %d\n",i);

			slicemultiply(&F, boxedbases[i], i, c);
			slicemultiply(&R, bases[i], i, c);
		}

		/* Now flatten R into a matrix */

		if (verbose)
			printf("\tFlattening residuals matrix...\n");

		Rmat = flatten_ndarray_to_sparse(&F, sidelen, 1, c);
				

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
				F.i[i] = malloc(data->rows * sizeof(int));
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
			int **oldi;

			oldi = F.i;
			F.i = malloc(F.ndim * sizeof(int *));
			for (i = 0; i < F.ndim; i++) {
				if (i % 2 == 0)
					F.i[i/2] = oldi[i];
				else
					F.i[F.ndim/2 + i/2] = oldi[i];
			}
			free(oldi);
		}

		/* Now flatten F */

		Fmat = flatten_ndarray_to_sparse(&F, sidelen, sidelen, c);
	
		/* XXX: optimization possibilities ended */

		scale2[0] = 1.0;
		fitmat = cholmod_add(Fmat, penalty, scale1, scale2, 1, 0, c);
		cholmod_free_sparse(&Fmat, c); /* we don't need Fmat anymore */

		Rdens = cholmod_sparse_to_dense(Rmat, c);
		cholmod_free_sparse(&Rmat, c); /* we don't need Fmat anymore */

		/*
		 * Now, we can solve the linear system 
		 */

		if (verbose)
			printf("Computing iteration %d least square solution...\n",n+1);
		
		coefficients = SuiteSparseQR_C_backslash_default(Fmat, Rdens, c);
	}
	
	/* Clean up detritus */

	if (verbose)
		printf("Done: cleaning up\n");

	for (i = 0; i < F.ndim; i++)
		free(F.i[i]);
	free(F.x); free(F.i); free(F.ranges);
	for (i = 0; i < R.ndim; i++)
		free(R.i[i]);
	free(R.x); free(R.i); free(R.ranges);

	cholmod_free_sparse(&Fmat, c);
	cholmod_free_dense(&Rdens, c);
	cholmod_free_sparse(&penalty, c);

	for (i = 0; i < out->ndim; i++) {
		cholmod_free_sparse(&bases[i], c);
		cholmod_free_sparse(&boxedbases[i], c);
	}
	free(bases); free(boxedbases);
	free(moduli);
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
	cholmod_free_dense(&coefficients, c);
}

static cholmod_sparse *
flatten_ndarray_to_sparse(struct ndsparse *array, size_t nrow, size_t ncol,
    cholmod_common *c)
{
	cholmod_triplet *trip;
	cholmod_sparse *sparse;
	int moduli[array->ndim];
	int i, j, k;

	trip = cholmod_allocate_triplet(nrow, ncol, array->rows, 0,
	    CHOLMOD_REAL, c);

	moduli[array->ndim-1] = 1;
	for (i = array->ndim-2; i >= 0; i--)
		moduli[i] = moduli[i+1]*array->ranges[i+1];

	for (i = 0; i < array->rows; i++) {
		k = 0;
		for (j = 0; j < array->ndim; j++)
			k += array->i[j][i]*moduli[j];

		((int *)(trip->j))[i] = k % ncol;
		((int *)(trip->i))[i] = k / nrow;
		((double *)(trip->x))[i] = array->x[i];
	}
	trip->nnz = array->rows;

	sparse = cholmod_triplet_to_sparse(trip, trip->nnz, c);
	cholmod_free_triplet(&trip, c);

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

	trip = cholmod_allocate_triplet(rows, cols,
	    max(rows,cols) - abs(diag), 0, CHOLMOD_REAL, c);

	while (row < rows && col < cols) {
		((int *)(trip->i))[trip->nnz] = row;
		((int *)(trip->j))[trip->nnz] = col;
		((double *)(trip->x))[trip->nnz] = 1.0;

		trip->nnz++; row++; col++;
	}

	sparse = cholmod_triplet_to_sparse(trip, trip->nnz, c);
	cholmod_free_triplet(&trip, c);

	return sparse;
}

cholmod_sparse *
calc_penalty(int *nsplines, int ndim, int dim, cholmod_common *c)
{
	cholmod_sparse *finitediff, *fd_trans, *DtD;
	cholmod_sparse *tmp, *tmp2, *result;
	double scale1[2] = {1.0, 0.0}, scale2[2] = {1.0, 0.0};
	int i;

	/* First, we will compute the finite difference matrix,
	 * which looks like this:
	 * 
	 * 1 -2  1  0 0 ...
	 * 0  1 -2  1 0 ...
	 * 0  0  1 -2 1 ...
	 */

	tmp = cholmod_speye(nsplines[dim] - 2, nsplines[dim],
	    CHOLMOD_REAL, c);
	tmp2 = offdiagonal(nsplines[dim] - 2, nsplines[dim],
	    1, c);
	scale2[0] = -2.0;
	finitediff = cholmod_add(tmp, tmp2, scale1, scale2, 1,
	    0, c);
	cholmod_free_sparse(&tmp2, c);
	cholmod_free_sparse(&tmp, c);

	scale2[0] = 1.0;
	tmp2 = offdiagonal(nsplines[dim] - 2, nsplines[dim],
	    2, c);
	tmp = cholmod_add(finitediff, tmp2, scale1, scale2, 1,
	    0, c);
	cholmod_free_sparse(&tmp2, c);
	cholmod_free_sparse(&finitediff, c);
	finitediff = tmp;

	/*
	 * Now we want DtD, which is the transpose of finitediff
	 * multiplied by finitediff
	 */

	fd_trans = cholmod_transpose(finitediff, 1, c);
	DtD = cholmod_ssmult(fd_trans, finitediff,
	    1 /* DtD is symmetric */, 1, 0, c);
	cholmod_free_sparse(&finitediff, c);
	cholmod_free_sparse(&fd_trans, c);


	/* Next take kronecker products to form the full P */

	tmp = NULL;
	result = NULL;
	for (i = 0; i < ndim; i++) {
		tmp2 = (i == dim) ? DtD : cholmod_speye(
		    nsplines[i], nsplines[i],
		    CHOLMOD_REAL, c);
		tmp2->stype = 1; /* The identity matrix is always symmetric. */

		if (result == NULL) {
			result = tmp2;
			continue;
		}

		tmp = kronecker_product(result, tmp2, c);
		cholmod_free_sparse(&result, c);
		cholmod_free_sparse(&tmp2, c);

		result = tmp;
	}

	return result;
}

