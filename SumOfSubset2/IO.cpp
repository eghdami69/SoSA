#define _CRT_SECURE_NO_DEPRECATE
#include "IO.h"
#include "mem.h"
#include"string.h"


/*
*	Read a 2d double matrix from file
*		matrix(O): output matrix
*		nr(I): number of rows
*		nc(I): number of columns
*		ifp(I): input file pointer
*/
void read_from_file_2d_double(double ** matrix, int nr, int nc, FILE * ifp){

	float temp;

	for (int i = 0; i < nr; i++){

		for (int j = 0; j < nc; j++){

			fscanf(ifp, "%f", &temp);
			matrix[i][j] = (double)temp;
		}

	}
}

/*
*	Read a double vector from file
*		matrix(O): output matrix
*		nc(I): number of columns
*		ifp(I): input file pointer
*/
void read_from_file_1d_double(double * matrix, int nc, FILE * ifp){

	float temp;

	for (int i = 0; i < nc; i++){

		fscanf(ifp, "%f", &temp);
		matrix[i] = (double)temp;

	}
}

void read_string_from_file(char * str, int n, FILE * ifp)
{

	if (fgets(str, n, ifp) == NULL) {

		str[0] = 'e';
	}

	int length = (int)strlen(str);

	// do this to eliminate the newline charachter ('\n') returned by fgets
	if ((str[0] != 'e') && (str[length - 1] == '\n')) {

		str[length - 1] = '\0';
	}

}

void read_from_file_2d_int(int ** array1, int nr, int nc, FILE * ofp){

	int j, k;
	int temp;

	for (j = 0; j < nr; j++){

		for (k = 0; k < nc; k++){

			fscanf(ofp, "%d", &temp);
			array1[j][k] = temp;
		}

	}
}
