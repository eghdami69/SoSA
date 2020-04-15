#ifndef MEM_H_
#define MEM_H_
#include <stdlib.h>

/* Allocate a double matrix */
void * allocmat(int rows, int columns, int size);
//void * f_matrix_calloc(int ndim, int size);
void * allocvec(int columns, int size);
void freemat(int **rp);
void freemat(double **rp);
//void freealn(alignment * aln);
void allocate_2d_int(int **& matrix, int nr, int nc);
void allocate_2d_double(double **& matrix, int nr, int nc);
void deallocate_2d_double(double **& matrix, int nr);
void deallocate_2d_int(int **& matrix, int nr);
void allocate_3d_int(int ***& matrix, int nx, int ny, int nz);
void deallocate_3d_int(int ***& matrix, int nx, int ny);

#endif /* MEM_H_ */
