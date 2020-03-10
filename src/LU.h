
#ifndef LU_H
#define LU_H

/******************************************************************************************/

#include <complex.h>
#include <math.h>
#include "ALLOCATION.h"

/******************************************************************************************/

//A small number.
#define TINY 1.0e-20;

/******************************************************************************************/

extern int ludcmp_double(double **a, int n, int *indx);

extern void lubksb_double(double **a, int n, int *indx, double *b);

extern void improve_double(double **a, double **alud, int n, int *indx, double *b, double *x);

extern int ludcmp_float(float **a, int n, int *indx);

extern void lubksb_float(float **a, int n, int *indx, float *b);

extern void improve_float(float **a, float **alud, int n, int *indx, float *b, float *x);

extern int ludcmp_complex(complex double **a, int n, int *indx);

extern void lubksb_complex(complex double **a, int n, int *indx, complex double *b);

extern void improve_complex(complex double **a, complex double **alud, int n, int *indx, \
                            complex double *b, complex double *x);

/******************************************************************************************/

#endif /* LU_H */
