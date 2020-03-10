
#ifndef ALLOCATION_H

#define ALLOCATION_H

/******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <complex.h>
#include <string.h>

/******************************************************************************************/

extern void nrerror(char error_text[]);

extern int *VECTOR_INT(long nl, long nh, int Flag_Initialize);

extern void FREE_VECTOR_INT(int *v, long nl);

extern float *VECTOR_FLOAT(long nl, long nh, int Flag_Initialize);

extern void FREE_VECTOR_FLOAT(float *v, long nl);

extern double *VECTOR_DOUBLE(long nl, long nh, int Flag_Initialize);

extern void FREE_VECTOR_DOUBLE(double *v, long nl);

extern complex double *VECTOR_COMPLEX(long nl, long nh, int Flag_Initialize);

extern void FREE_VECTOR_COMPLEX(complex double *v, long nl);

extern char **MATRIX_CHAR(long nrl, long nrh, long ncl, long nch, int Flag_Initialize);

extern void FREE_MATRIX_CHAR(char **m, long nrl, long ncl);

extern int **MATRIX_INT(long nrl, long nrh, long ncl, long nch, int Flag_Initialize);

extern void FREE_MATRIX_INT(int **m, long nrl, long ncl);

extern float **MATRIX_FLOAT(long nrl, long nrh, long ncl, long nch, int Flag_Initialize);

extern void FREE_MATRIX_FLOAT(float **m, long nrl, long ncl);

extern double **MATRIX_DOUBLE(long nrl, long nrh, long ncl, long nch, int Flag_Initialize);

extern void FREE_MATRIX_DOUBLE(double **m, long nrl, long ncl);

extern complex double **MATRIX_COMPLEX(long nrl, long nrh, long ncl, long nch, \
                                       int Flag_Initialize);

extern void FREE_MATRIX_COMPLEX(complex double **m, long nrl, long ncl);

extern float **MATRIX_DIAGONAL(long nh, int Flag_Initialize);

extern void FREE_MATRIX_DIAGONAL(float **m);

extern complex double **MATRIX_RHO(long nh, int Flag_Initialize);

extern void FREE_MATRIX_RHO(complex double **Rho);

extern complex double ***MATRIX3_RHO(long n1, long n2, int Flag_Intialize);

extern void FREE_MATRIX3_RHO(complex double ***Rho);

extern double ***MATRIX3_RHO_DB(long n1, long n2, int Flag_Initialize);

extern void FREE_MATRIX3_RHO_DB(double ***Rho);

extern double ***CUBE_DOUBLE(long Ni0, long Ni1, long Nj0, long Nj1, long Nk0, long Nk1, int Flag_Initialize);

extern void FREE_CUBE_DOUBLE(double ***C, long Ni0, long Nj0, long Nk0);

/******************************************************************************************/

#endif /* ALLOCATION_H */
