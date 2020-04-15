#include "IO.h"
#include "mem.h"



/* Allocate vector */
void *allocvec(int columns, int size)
{
	void *p;
	p = calloc(columns, size);
	return p;
}



//allocate matrix
void * allocmat(int rows, int columns, int size)
{
	int             i;
	void          **rp;

	rp = (void **)malloc(rows * sizeof(void *));

	for (i = 0; i < rows; i++)
		rp[i] = calloc(columns, size);

	return rp;
}

//free matrix
void freemat(int **rp)
{
	int rows;
	rows = sizeof(rp) / sizeof(*rp);
	while (rows--)
		free(rp[rows]);
	free(rp);
}

void freemat(double **rp)
{
	int rows;
	rows = sizeof(rp) / sizeof(*rp);
	while (rows--)
		free(rp[rows]);
	free(rp);
}

/*
*	Allocate space to a 2d double matrix
*		matrix(O): output matrix
*		nr(I): number of rows
*		nc(I): number of columns
*/
void allocate_2d_double(double **& matrix, int nr, int nc)
{
	matrix = new double*[nr];

	for (int i = 0; i < nr; i++) {

		matrix[i] = new double[nc];

		for (int j = 0; j < nc; j++) {

			matrix[i][j] = 0.0;
		}
	}
}

/*
*	Allocate space to a 2d int matrix
*		matrix(O): output matrix
*		nr(I): number of rows
*		nc(I): number of columns
*/
void allocate_2d_int(int **& matrix, int nr, int nc)
{
	matrix = new int*[nr];

	for (int i = 0; i < nr; i++) {

		matrix[i] = new int[nc];

		for (int j = 0; j < nc; j++) {

			matrix[i][j] = 0;
		}
	}
}

void deallocate_2d_double(double **& matrix, int nr)
{

	for (int j = 0; j < nr; j++)

		delete[] matrix[j];

	delete[] matrix;

}

/*
*	Deallocate space of a 2d int matrix
*		matrix(I): input matrix
*		nr(I): number of rows
*/
void deallocate_2d_int(int **& matrix, int nr)
{

	for (int j = 0; j < nr; j++)

		delete[] matrix[j];

	delete[] matrix;

}

void allocate_3d_int(int ***& matrix, int nx, int ny, int nz)
{

	matrix = new int**[nx];

	for (int i = 0; i < nx; i++) {

		matrix[i] = new int*[ny];

		for (int j = 0; j < ny; j++) {

			matrix[i][j] = new int[nz];

			for (int k = 0; k < nz; k++) {

				matrix[i][j][k] = 0;
			}
		}
	}

}

void deallocate_3d_int(int ***& matrix, int nx, int ny)
{

	int j, k;

	for (j = 0; j < nx; j++) {

		for (k = 0; k < ny; k++)

			delete[] matrix[j][k];

	}

	for (j = 0; j < nx; j++)

		delete[] matrix[j];

	delete[] matrix;

}
