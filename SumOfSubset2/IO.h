#ifndef IO2_H_
#define IO2_H_
#include <stdio.h>

void read_from_file_2d_double(double ** matrix, int nr, int nc, FILE * ifp);
void read_from_file_1d_double(double * matrix, int nc, FILE * ifp);
void read_string_from_file(char * str, int n, FILE * ifp);
void read_from_file_2d_int(int ** array1, int nr, int nc, FILE * ofp);

#endif /* IO_H_ */
