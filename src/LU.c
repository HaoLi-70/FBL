
#include "LU.h"

/**********************************************************************************************/

extern int ludcmp_double(double **a, int n, int *indx){
    
    /******************************************************************************************
     Purpose:
     LU decomposition (double matrix).
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     a[1..n][1..n], the input matrix.
     n, the size of the matrix.
     Output parameters:
     a[1..n][1..n], the output matrix.
     indx[1..n], records the row permutation.
     Return:
     return ±1 depending on whether the number of row interchanges was even or odd, respectively.
     Reference:
     Numerical recipes in C 2ed.
     Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
     permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14)
     above; indx[1..n] is an output vector that records the row permutation effected by the
     partial pivoting; d is output as ±1 depending on whether the number of row interchanges
     was even or odd, respectively. This routine is used in combination with lubksb to solve
     linear equations or invert a matrix.
     ******************************************************************************************/
    
    int i, j, k, d, imax = 0;
    double big,dum,sum,temp;
    double *vv;
    vv = VECTOR_DOUBLE(1, n, 1);
    d = 1;
    
    for (i=1;i<=n;i++) {
        big = 0.0;
        for (j=1;j<=n;j++)
            if ((temp = fabs(a[i][j])) > big) big=temp;
        if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
        vv[i] = 1.0/big;
    }
    
    for (j=1;j<=n;j++) {
        for (i=1;i<j;i++) {
            sum = a[i][j];
            for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
        }
        big = 0.0;
        for (i=j;i<=n;i++) {
            sum = a[i][j];
            for (k=1;k<j;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
            if ( (dum=vv[i]*fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (k=1;k<=n;k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0) a[j][j] = TINY;
        if (j != n) {
            dum = 1.0/(a[j][j]);
            for (i=j+1;i<=n;i++) a[i][j] *= dum;
        }
    }
    
    FREE_VECTOR_DOUBLE(vv,1);
    return d;
}


void lubksb_double(double **a, int n, int *indx, double *b){
    
    /******************************************************************************************
     Purpose:
     Solve linear equations (double matrix).
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     a[1..n][1..n], the input matrix (LU decomposition).
     n, the size of the matrix.
     indx[1..n], records the row permutation.
     b[1..n], the right-hand side vector.
     Output parameters:
     b[1..n], the solutions.
     Reference:
     Numerical recipes in C 2ed.
     Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the matrix
     A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
     as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
     B, and returns with the solution vector X. a, n, and indx are not modified by this routine
     and can be left in place for successive calls with different right-hand sides b. This routine
     takes into account the possibility that b will begin with many zero elements, so it is
     efficient for use in matrix inversion.
     ******************************************************************************************/

    int i, ip, j, ii = 0;
    double sum;
    
    for (i=1;i<=n;i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii){
            for (j=ii;j<=i-1;j++)
                sum -= a[i][j]*b[j];
        }
        else if (sum) {
            ii = i;
        }
        b[i] = sum;
    }
    
    for (i=n;i>=1;i--) {
        sum = b[i];
        for (j=i+1;j<=n;j++){
            sum -= a[i][j]*b[j];
        }
        b[i] = sum/a[i][i];
    }
    return;
}

extern void improve_double(double **a, double **alud, int n, int *indx, double *b, double *x){
    
    /******************************************************************************************
     Purpose:
     Improves a solution vector x[1..n] of the linear set of equations A · X = B (double matrix).
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     a[1..n][1..n], the A matrix.
     alud[1..n][1..n], the LU decomposition matrix.
     n, the size of the matrix.
     indx[1..n], records the row permutation.
     b[1..n], the right-hand side vector.
     x[1..n], the solutions.
     Output parameters:
     x[1..n], the imporved solutions.
     Reference:
     Numerical recipes in C 2ed.
     Improves a solution vector x[1..n] of the linear set of equations A · X = B. The matrix
     a[1..n][1..n], and the vectors b[1..n] and x[1..n] are input, as is the dimension n. Also
     input is alud[1..n][1..n], the LU decomposition of a as returned by ludcmp, and the vector
     indx[1..n] also returned by that routine. On output, only x[1..n] is modified, to an
     improved set of values.
     ******************************************************************************************/
    
    int i, j;
    double sdp, *r;
    r = VECTOR_DOUBLE(1,n,1);
    
    for (i=1;i<=n;i++) {
        sdp = -b[i];
        for (j=1;j<=n;j++){
            sdp += a[i][j]*x[j];
        }
        r[i] = sdp;
    }
    
    lubksb_double(alud,n,indx,r);
    
    for (i=1;i<=n;i++)
        x[i] -= r[i];
    
    FREE_VECTOR_DOUBLE(r,1);
    return;
}

extern int ludcmp_float(float **a, int n, int *indx){
    
    /******************************************************************************************
     Purpose:
     LU decomposition (float matrix).
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     a[1..n][1..n], the input matrix.
     n, the size of the matrix.
     Output parameters:
     a[1..n][1..n], the output matrix.
     indx[1..n], records the row permutation.
     Return:
     return ±1 depending on whether the number of row interchanges was even or odd, respectively.
     Reference:
     Numerical recipes in C 2ed.
     Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
     permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14)
     above; indx[1..n] is an output vector that records the row permutation effected by the
     partial pivoting; d is output as ±1 depending on whether the number of row interchanges
     was even or odd, respectively. This routine is used in combination with lubksb to solve
     linear equations or invert a matrix.
     ******************************************************************************************/
    
    int i, j, k, d, imax = 0;
    float big, dum, sum, temp;
    double *vv;
    vv = VECTOR_DOUBLE(1,n,1);
    d = 1;
    
    for (i=1;i<=n;i++) {
        big = 0.0;
        for (j=1;j<=n;j++){
            if ((temp=fabs(a[i][j])) > big){
                big = temp;
            }
        }
        if (big == 0.0)
            nrerror("Singular matrix in routine ludcmp");
        vv[i] = 1.0/big;
    }
    
    for (j=1;j<=n;j++) {
        for (i=1;i<j;i++) {
            sum = a[i][j];
            for (k=1;k<i;k++){
                sum -= a[i][k]*a[k][j];
            }
            a[i][j] = sum;
        }
        big = 0.0;
        for (i=j;i<=n;i++) {
            sum = a[i][j];
            for (k=1;k<j;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
            if ( (dum=vv[i]*fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (k=1;k<=n;k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0)
            a[j][j] = TINY;
        if (j != n) {
            dum = 1.0/(a[j][j]);
            for (i=j+1;i<=n;i++)
                a[i][j] *= dum;
        }
    }
    
    FREE_VECTOR_DOUBLE(vv,1);
    return d;
}


void lubksb_float(float **a, int n, int *indx, float *b){
    
    /******************************************************************************************
     Purpose:
     Solve linear equations (float matrix).
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     a[1..n][1..n], the input matrix (LU decomposition).
     n, the size of the matrix.
     indx[1..n], records the row permutation.
     b[1..n], the right-hand side vector.
     Output parameters:
     b[1..n], the solutions.
     Reference:
     Numerical recipes in C 2ed.
     Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the matrix
     A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
     as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
     B, and returns with the solution vector X. a, n, and indx are not modified by this routine
     and can be left in place for successive calls with different right-hand sides b. This routine
     takes into account the possibility that b will begin with many zero elements, so it is
     efficient for use in matrix inversion.
     ******************************************************************************************/
    
    int i, ip, j, ii = 0;
    float sum;
    
    for (i=1;i<=n;i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii){
            for (j=ii;j<=i-1;j++){
                sum -= a[i][j]*b[j];
            }
        }else if (sum){
            ii = i;
        }
        b[i] = sum;
    }
    
    for (i=n;i>=1;i--) {
        sum = b[i];
        for (j=i+1;j<=n;j++){
            sum -= a[i][j]*b[j];
        }
        b[i] = sum/a[i][i];
    }
    
    return;
}

extern void improve_float(float **a, float **alud, int n, int *indx, float *b, float *x){
    
    /******************************************************************************************
     Purpose:
     Improves a solution vector x[1..n] of the linear set of equations A · X = B (float matrix).
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     a[1..n][1..n], the A matrix.
     alud[1..n][1..n], the LU decomposition matrix.
     n, the size of the matrix.
     indx[1..n], records the row permutation.
     b[1..n], the right-hand side vector.
     x[1..n], the solutions.
     Output parameters:
     x[1..n], the imporved solutions.
     Reference:
     Numerical recipes in C 2ed.
     Improves a solution vector x[1..n] of the linear set of equations A · X = B. The matrix
     a[1..n][1..n], and the vectors b[1..n] and x[1..n] are input, as is the dimension n. Also
     input is alud[1..n][1..n], the LU decomposition of a as returned by ludcmp, and the vector
     indx[1..n] also returned by that routine. On output, only x[1..n] is modified, to an
     improved set of values.
     ******************************************************************************************/
    
    int i, j;
    float sdp,  *r;
    r = VECTOR_FLOAT(1,n,1);
    
    for (i=1;i<=n;i++) {
        sdp = -b[i];
        for (j=1;j<=n;j++){
            sdp += a[i][j]*x[j];
        }
        r[i] = sdp;
    }
    
    lubksb_float(alud,n,indx,r);
    
    for (i=1;i<=n;i++)
        x[i] -= r[i];
    
    FREE_VECTOR_FLOAT(r,1);
    return;
}


extern int ludcmp_complex(complex double **a, int n, int *indx){
    
    /******************************************************************************************
     Purpose:
     LU decomposition (complex double matrix).
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     a[1..n][1..n], the input matrix.
     n, the size of the matrix.
     Output parameters:
     a[1..n][1..n], the output matrix.
     indx[1..n], records the row permutation.
     Return:
     return ±1 depending on whether the number of row interchanges was even or odd, respectively.
     Reference:
     Numerical recipes in C 2ed.
     Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
     permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14)
     above; indx[1..n] is an output vector that records the row permutation effected by the
     partial pivoting; d is output as ±1 depending on whether the number of row interchanges
     was even or odd, respectively. This routine is used in combination with lubksb to solve
     linear equations or invert a matrix.
     ******************************************************************************************/
    
    int i, j, k, d, imax = 0;
    double big, temp, *vv;
    complex double sum,dum;
    vv = VECTOR_DOUBLE(1,n,1);
    d = 1;
    
    for (i=1;i<=n;i++) {
        big = 0.0;
        for (j=1;j<=n;j++)
            if ((temp=creal(a[i][j])*creal(a[i][j])+cimag(a[i][j])*cimag(a[i][j])) > big)
                big = temp;
        if (big == 0.0)
            nrerror("Singular matrix in routine ludcmp");
        vv[i] = 1.0/big;
    }
    
    for (j=1;j<=n;j++) {
        for (i=1;i<j;i++) {
            sum = a[i][j];
            for (k=1;k<i;k++){
                sum -= a[i][k]*a[k][j];
            }
            a[i][j] = sum;
        }
        big = 0.0;
        for (i=j;i<=n;i++) {
            sum = a[i][j];
            for (k=1;k<j;k++){
                sum -= a[i][k]*a[k][j];
            }
            a[i][j] = sum;
            if ( (temp=vv[i]*(creal(sum)*creal(sum)+cimag(sum)*cimag(sum))) >= big) {
                big = temp;
                imax = i;
            }
        }
        if (j != imax) {
            for (k=1;k<=n;k++) {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0) a[j][j]=TINY;
        if (j != n) {
            dum = 1.0/(a[j][j]);
            for (i=j+1;i<=n;i++)
                a[i][j] *= dum;
        }
    }
    
    FREE_VECTOR_DOUBLE(vv,1);
    return d;
}


extern void lubksb_complex(complex double **a, int n, int *indx, complex double *b){
    
    /******************************************************************************************
     Purpose:
     Solve linear equations (complex double matrix).
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     a[1..n][1..n], the input matrix (LU decomposition).
     n, the size of the matrix.
     indx[1..n], records the row permutation.
     b[1..n], the right-hand side vector.
     Output parameters:
     b[1..n], the solutions.
     Reference:
     Numerical recipes in C 2ed.
     Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the matrix
     A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
     as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
     B, and returns with the solution vector X. a, n, and indx are not modified by this routine
     and can be left in place for successive calls with different right-hand sides b. This routine
     takes into account the possibility that b will begin with many zero elements, so it is
     efficient for use in matrix inversion.
     ******************************************************************************************/
    
    int i, ip, j, ii = 0;
    complex double sum;
    
    for (i=1;i<=n;i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii){
            for (j=ii;j<=i-1;j++)
                sum -= a[i][j]*b[j];
        }
        else if (sum){
            ii = i;
        }
        b[i] = sum;
    }
    for (i=n;i>=1;i--) {
        sum = b[i];
        for (j=i+1;j<=n;j++){
            sum -= a[i][j]*b[j];
        }
        b[i] = sum/a[i][i];
    }
    
    return;
}


extern void improve_complex(complex double **a, complex double **alud, int n, int *indx, \
                            complex double *b, complex double *x){
    
    /******************************************************************************************
     Purpose:
     Improves a solution vector x[1..n] of the linear set of equations A · X = B (complex matrix).
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     a[1..n][1..n], the A matrix.
     alud[1..n][1..n], the LU decomposition matrix.
     n, the size of the matrix.
     indx[1..n], records the row permutation.
     b[1..n], the right-hand side vector.
     x[1..n], the solutions.
     Output parameters:
     x[1..n], the imporved solutions.
     Reference:
     Numerical recipes in C 2ed.
     Improves a solution vector x[1..n] of the linear set of equations A · X = B. The matrix
     a[1..n][1..n], and the vectors b[1..n] and x[1..n] are input, as is the dimension n. Also
     input is alud[1..n][1..n], the LU decomposition of a as returned by ludcmp, and the vector
     indx[1..n] also returned by that routine. On output, only x[1..n] is modified, to an
     improved set of values.
     ******************************************************************************************/
    
    int i, j;
    complex double sdp, *r;
    r = VECTOR_COMPLEX(1,n,1);
    
    for (i=1;i<=n;i++) {
        sdp = -b[i];
        for (j=1;j<=n;j++){
            sdp += a[i][j]*x[j];
        }
        r[i] = sdp;
    }
    
    lubksb_complex(alud,n,indx,r);
    
    for (i=1;i<=n;i++)
        x[i] -= r[i];
    
    FREE_VECTOR_COMPLEX(r,1);
    return;
}

/**********************************************************************************************/

