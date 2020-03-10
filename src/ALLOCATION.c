
#include "ALLOCATION.h"

/**********************************************************************************************/

extern void nrerror(char error_text[]){
    
    /******************************************************************************************
     Purpose:
     Numerical Recipes standard error handler.
     Record of revisions:
     21 Mar 2018
     Input parameters:
     error_text[], the error text.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    fprintf(stderr,"%s\n",error_text);
    
    fprintf(stderr,"...now exiting to system...\n");
    
    exit(1);
    
    return;
}


extern int *VECTOR_INT(long nl, long nh, int Flag_Initialize){
    
    /******************************************************************************************
     Purpose:
     Allocate an int vector with subscript range v[nl..nh].
     Record of revisions:
     7 Jun. 2019
     Input parameters:
     nl, the first index of subscript range.
     nh, the last index of subscript range.
     Return:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    int *v;

    long Byte = (nh-nl+1)*sizeof(int);
    
    v = (int *)malloc(Byte);
    
    if(Flag_Initialize)
        memset(v,0,Byte);
  
    if (!v)
        nrerror("allocation failure in VECTOR_INT()");
    
    return v-nl;
}


extern void FREE_VECTOR_INT(int *v, long nl){
    
    /******************************************************************************************
     Purpose:
     Free an int vector allocated with VECTOR_INT().
     Record of revisions:
     21 Mar 2018
     Input parameters:
     *v, the pointer.
     nl, the first index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    free((v+nl));
    
    return;
}


extern float *VECTOR_FLOAT(long nl, long nh, int Flag_Initialize){
    
    /******************************************************************************************
     Purpose:
     Allocate a float vector with subscript range v[nl..nh].
     Record of revisions:
     7 Jun. 2019
     Input parameters:
     nl, the first index of subscript range.
     nh, the last index of subscript range.
     Return:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    float *v;
    
    long Byte = (nh-nl+1)*sizeof(float);
    
    v = (float *)malloc(Byte);
    
    if(Flag_Initialize)
        memset(v, 0, Byte);

    if (!v)
        nrerror("allocation failure in VECTOR_FLOAT()");
    
    return v-nl;
}


extern void FREE_VECTOR_FLOAT(float *v, long nl){
    
    /******************************************************************************************
     Purpose:
     Free a float vector allocated with VECTOR_FLOAT().
     Record of revisions:
     21 Mar 2018
     Input parameters:
     *v, the pointer.
     nl, the first index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    free((v+nl));
    
    return;
}


extern double *VECTOR_DOUBLE(long nl, long nh, int Flag_Initialize){
    
    /******************************************************************************************
     Purpose:
     Allocate a double vector with subscript range v[nl..nh].
     Record of revisions:
     7 Jun. 2019
     Input parameters:
     nl, the first index of subscript range.
     nh, the last index of subscript range.
     Return:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    double *v;
    
    long Byte=(nh-nl+1)*sizeof(double);
    
    v=(double *)malloc(Byte);
    
    if (Flag_Initialize)
        memset(v, 0, Byte);

    if (!v)
        nrerror("allocation failure in VECTOR_DOUBLE()");
    
    return v-nl;
}


extern void FREE_VECTOR_DOUBLE(double *v, long nl){
    
    /******************************************************************************************
     Purpose:
     Free a double vector allocated with VECTOR_DOUBLE().
     Record of revisions:
     21 Mar 2018
     Input parameters:
     *v, the pointer.
     nl, the first index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    free((v+nl));
    
    return;
}


extern complex double *VECTOR_COMPLEX(long nl, long nh, int Flag_Initialize){
    
    /******************************************************************************************
     Purpose:
     Allocate a complex double vector with subscript range v[nl..nh].
     Record of revisions:
     7th Jun. 2019
     Input parameters:
     nl, the first index of subscript range.
     nh, the last index of subscript range.
     Return:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    complex double *v;
    
    long Byte = (nh-nl+1)*sizeof(complex double);
    
    v = (complex double *)malloc(Byte);
    
    if (Flag_Initialize)
        memset(v, 0, Byte);
    
    if (!v)
        nrerror("allocation failure in VECTOR_COMPLEX()");
    
    return v-nl;
}


extern void FREE_VECTOR_COMPLEX(complex double *v, long nl){
    
    /******************************************************************************************
     Purpose:
     Free a complex double vector allocated with VECTOR_COMPLEX().
     Record of revisions:
     21 Mar 2018
     Input parameters:
     *v, the pointer.
     nl, the first index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    free((v+nl));
    
    return;
}


extern char **MATRIX_CHAR(long nrl, long nrh, long ncl, long nch, int Flag_Initialize){
    
    /******************************************************************************************
     Purpose:
     Allocate an char matrix with subscript range m[nrl..nrh][ncl..nch].
     Record of revisions:
     7th Jun. 2019
     Input parameters:
     nrl, the first row index of subscript range.
     nrh, the last row index of subscript range.
     ncl, the first column index of subscript range.
     nch, the last column of subscript range.
     Return:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
    
    char **m;
    
    m = (char **) malloc(((nrow)*sizeof(char*)));
    
    if (!m)
        nrerror("allocation failure 1 in MATRIX_CHAR()");
    
    m -= nrl;
    
    long Byte = (nrow*ncol)*sizeof(char);
    
    m[nrl] = (char *) malloc(Byte);
    
    if (Flag_Initialize)
        memset(m[nrl], 0, Byte);
    
    if (!m[nrl])
        
        nrerror("allocation failure 2 in MATRIX_CHAR()");
    
    m[nrl] -= ncl;
    
    for(i=nrl+1; i<=nrh; i++)
        
        m[i] = m[i-1]+ncol;
    
    return m;
}


extern void FREE_MATRIX_CHAR(char **m, long nrl, long ncl){
    
    /******************************************************************************************
     Purpose:
     Free a char matrix allocated with MATRIX_CHAR().
     Record of revisions:
     21 Mar 2018
     Input parameters:
     **m, the pointer.
     nrl, the first row index of subscript range.
     ncl, the first column index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    free(m[nrl]+ncl);
    
    free(m+nrl);
    
    return;
}


extern int **MATRIX_INT(long nrl, long nrh, long ncl, long nch, int Flag_Initialize){
    
    /******************************************************************************************
     Purpose:
     Allocate an int matrix with subscript range m[nrl..nrh][ncl..nch].
     Record of revisions:
     7th Jun. 2019
     Input parameters:
     nrl, the first row index of subscript range.
     nrh, the last row index of subscript range.
     ncl, the first column index of subscript range.
     nch, the last column of subscript range.
     Return:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
    
    int **m;
    
    m = (int **) malloc(((nrow)*sizeof(int*)));
    
    if (!m)
        nrerror("allocation failure 1 in MATRIX_INT()");
    
    m -= nrl;
    
    long Byte = (nrow*ncol)*sizeof(int);
    
    m[nrl] = (int *) malloc(Byte);
    
    if (Flag_Initialize)
        memset(m[nrl], 0, Byte);
    
    if (!m[nrl])
        nrerror("allocation failure 2 in MATRIX_INT()");
    
    m[nrl] -= ncl;
    
    for(i=nrl+1; i<=nrh; i++)
        m[i] = m[i-1]+ncol;
    
    return m;
}


extern void FREE_MATRIX_INT(int **m, long nrl, long ncl){
    
    /******************************************************************************************
     Purpose:
     Free an int matrix allocated with MATRIX_INT().
     Record of revisions:
     21 Mar 2018
     Input parameters:
     **m, the pointer.
     nrl, the first row index of subscript range.
     ncl, the first column index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    free(m[nrl]+ncl);
    
    free(m+nrl);
    
    return;
}


extern float **MATRIX_FLOAT(long nrl, long nrh, long ncl, long nch, int Flag_Initialize){
    
    /******************************************************************************************
     Purpose:
     Allocate a float matrix with subscript range m[nrl..nrh][ncl..nch].
     Record of revisions:
     7th Jun 2019
     Input parameters:
     nrl, the first row index of subscript range.
     nrh, the last row index of subscript range.
     ncl, the first column index of subscript range.
     nch, the last column of subscript range.
     Return:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
    
    float **m;
    
    m = (float **) malloc(((nrow)*sizeof(float*)));
    
    if (!m)
        nrerror("allocation failure 1 in MATRIX_FLOAT()");
    
    m -= nrl;

    long Byte = (nrow*ncol)*sizeof(float);
    
    m[nrl] = (float *) malloc(Byte);
    
    if (Flag_Initialize)
        memset(m[nrl], 0, Byte);
   
    if (!m[nrl])
        nrerror("allocation failure 2 in MATRIX_FLOAT()");
    
    m[nrl] -= ncl;
    
    for(i=nrl+1; i<=nrh; i++)
        m[i] = m[i-1]+ncol;
    
    return m;
}


extern void FREE_MATRIX_FLOAT(float **m, long nrl, long ncl){
    
    /******************************************************************************************
     Purpose:
     Free a float matrix allocated with MATRIX_FLOAT().
     Record of revisions:
     21 Mar 2018
     Input parameters:
     **m, the pointer.
     nrl, the first row index of subscript range.
     ncl, the first column index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    free(m[nrl]+ncl);
    
    free(m+nrl);
    
    return;
}


extern double **MATRIX_DOUBLE(long nrl, long nrh, long ncl, long nch, int Flag_Initialize){
    
    /******************************************************************************************
     Purpose:
     Allocate a double matrix with subscript range m[nrl..nrh][ncl..nch].
     Record of revisions:
     7rh Jun. 2019
     Input parameters:
     nrl, the first row index of subscript range.
     nrh, the last row index of subscript range.
     ncl, the first column index of subscript range.
     nch, the last column of subscript range.
     Return:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
    
    double **m;
    
    m = (double **) malloc(((nrow)*sizeof(double*)));
    
    if (!m)
        nrerror("allocation failure 1 in MATRIX_DOUBLE()");
    
    m -= nrl;

    long Byte = (nrow*ncol)*sizeof(double);
    
    m[nrl] = (double *) malloc(Byte);
    
    if (Flag_Initialize)
        memset(m[nrl], 0, Byte);
    
    if (!m[nrl])
        nrerror("allocation failure 2 in MATRIX_DOUBLE()");
    
    m[nrl] -= ncl;
    
    for(i=nrl+1; i<=nrh; i++)
        m[i] = m[i-1]+ncol;
    
    return m;
}


extern void FREE_MATRIX_DOUBLE(double **m, long nrl, long ncl){
    
    /******************************************************************************************
     Purpose:
     Free a double matrix allocated with MATRIX_DOUBLE().
     Record of revisions:
     21 Mar 2018
     Input parameters:
     **m, the pointer.
     nrl, the first row index of subscript range.
     ncl, the first column index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    free((m[nrl]+ncl));
    
    free((m+nrl));
    
    return;
}


extern complex double **MATRIX_COMPLEX(long nrl, long nrh, long ncl, long nch, int Flag_Initialize){
    
    /******************************************************************************************
     Purpose:
     Allocate a complex double matrix with subscript range m[nrl..nrh][ncl..nch].
     Record of revisions:
     7th Jun. 2019
     Input parameters:
     nrl, the first row index of subscript range.
     nrh, the last row index of subscript range.
     ncl, the first column index of subscript range.
     nch, the last column of subscript range.
     Return:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
    
    complex double **m;
    
    m = (complex double **) malloc(((nrow)*sizeof(complex double*)));
    
    if (!m)
        nrerror("allocation failure 1 in MATRIX_COMPLEX()");
    
    m -= nrl;
    
    long Byte = (nrow*ncol)*sizeof(complex double);
    
    m[nrl] = (complex double *) malloc(Byte);
    
    if (Flag_Initialize)
        memset(m[nrl], 0, Byte);
   
    if (!m[nrl])
        nrerror("allocation failure 2 in MATRIX_COMPLEX()");
    
    m[nrl] -= ncl;
    
    for(i=nrl+1; i<=nrh; i++)
        m[i] = m[i-1]+ncol;
    
    return m;
}


extern void FREE_MATRIX_COMPLEX(complex double **m, long nrl, long ncl){
    
    /******************************************************************************************
     Purpose:
     Free a complex double matrix allocated with MATRIX_COMPLEX().
     Record of revisions:
     21 Mar 2018
     Input parameters:
     **m, the pointer.
     nrl, the first row index of subscript range.
     ncl, the first column index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    free((m[nrl]+ncl));
    
    free((m+nrl));
    
    return;
}


extern float **MATRIX_DIAGONAL(long nh, int Flag_Initialize){
    
    /******************************************************************************************
     Purpose:
     Allocate a float diagonal matrix with subscript range m[0..nh][0..(nh)].
     Record of revisions:
     7th Jun. 2019
     Input parameters:
     nh, the dimension of the diagonal matrix..
     Return:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    long nrow = nh+1,i;
    
    float **m;
    
    m = (float **) malloc(nrow*sizeof(float*));
    
    if (!m)
        nrerror("allocation failure 1 in MATRIX_DIAGONAL()");
    
    long Byte = (nh+2)*nrow/2*sizeof(float);
    
    m[0] = (float *) malloc(Byte);
    
    if (Flag_Initialize)
        memset(m[0], 0, Byte);
    
    if (!m[0])
        nrerror("allocation failure 2 in MATRIX_DIAGONAL()");
    
    for(i=1; i<=nh; i++)
        m[i] = m[i-1]+i;
    
    return m;

}


extern void FREE_MATRIX_DIAGONAL(float **m){
    
    /******************************************************************************************
     Purpose:
     Free an float diagonal matrix allocated with MATRIX_DIAGONAL().
     Record of revisions:
     21 Mar 2018
     Input parameters:
     **m, the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    free(m[0]);
    
    free(m);
    
    return;
}


extern complex double **MATRIX_RHO(long nh, int Flag_Initialize){
    
    /******************************************************************************************
     Purpose:
     Allocate an complex double rho matrix with subscript range m[0..nh][0..(nh*2+1)].
     Record of revisions:
     7th Jun. 2019
     Input parameters:
     nh, the dimension of the rho matrix.
     Return:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    long nrow = nh+1,i;
    
    complex double **R;
    
    R = (complex double **)malloc(nrow*sizeof(complex double *));
    
    if (!R)
        nrerror("allocation failure 1 in MATRIX_RHO()");
    
    long Byte = nrow*nrow*sizeof(complex double);
    
    R[0] = (complex double *)malloc(Byte);
    
    if (Flag_Initialize)
        memset(R[0], 0, Byte);
   
    if (!R[0])
        nrerror("allocation failure 2 in MATRIX_RHO()");
    
    for (i=1; i<=nh; i++)
        R[i] = R[i-1]+2*i;
    
    return R;
}


extern void FREE_MATRIX_RHO(complex double **Rho){
    
    /******************************************************************************************
     Purpose:
     Free a complex rho matrix allocated by MATRIX_RHO().
     Record of revisions:
     21 Mar 2018
     Input parameters:
     **Rho, the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    free((Rho[0]));
    
    free((Rho));
    
    return;
    
}


extern complex double ***MATRIX3_RHO(long n1, long n2, int Flag_Initialize){
    
    /******************************************************************************************
     Purpose:
     Allocate an complex double rho matrix m[i][j][k] with subscript range 0<=i<=n1,
     0<=j<=n2, -j<=k<=j.
     Record of revisions:
     7th Jun. 2019
     Input parameters:
     nh, the dimension of the rho matrix.
     Return:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    long nx = n1+1, ny = n2+1, i, j;
    
    complex double ***R;
    
    R = (complex double ***)malloc(nx*sizeof(complex double **));
    
    if (!R)
        nrerror("allocation failure 1 in MATRIX3_RHO()");
    
    long Bt1 = nx*ny*sizeof(complex double *);
    
    R[0] = (complex double **)malloc(Bt1);
    
    if (!R[0])
        nrerror("allocation failure 2 in MATRIX3_RHO()");
    
    for (i=1; i<nx; i++)
        R[i] = R[i-1]+ny;
    
    long Bt2 = nx*ny*ny*sizeof(complex double);
    
    R[0][0] = (complex double *)malloc(Bt2);
    
    if (!R[0][0])
        nrerror("allocation failure 3 in MATRIX3_RHO()");
    
    if (Flag_Initialize)
        memset(R[0][0], 0, Bt2);
    
    for (i=1; i<nx; i++) {
        R[i][0] = R[i-1][0]+ny*ny;
    }
    
    for (i=0; i<nx; i++) {
        for (j=1; j <ny; j++) {
            R[i][j] = R[i][j-1]+2*j;
        }
    }
    
    return R;
}


extern void FREE_MATRIX3_RHO(complex double ***Rho){
    
    /******************************************************************************************
     Purpose:
     Free a complex rho matrix allocated by MATRIX3_RHO().
     Record of revisions:
     7th Jun. 2019
     Input parameters:
     ***Rho, the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    free(Rho[0][0]);
    
    free(Rho[0]);
    
    free(Rho);
    
    return;
}

extern double ***MATRIX3_RHO_DB(long n1, long n2, int Flag_Initialize){
    
    /******************************************************************************************
     Purpose:
     Allocate an complex double rho matrix m[i][j][k] with subscript range 0<=i<=n1,
     0<=j<=n2, -j<=k<=j.
     Record of revisions:
     7th Jun. 2019
     Input parameters:
     nh, the dimension of the rho matrix.
     Return:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    long nx = n1+1, ny = n2+1, i, j;
    
    double ***R;
    
    R = (double ***)malloc(nx*sizeof(double **));
    
    if (!R)
        nrerror("allocation failure 1 in MATRIX3_RHO()");
    
    long Bt1 = nx*ny*sizeof(double *);
    
    R[0] = (double **)malloc(Bt1);
    
    if (!R[0])
        nrerror("allocation failure 2 in MATRIX3_RHO()");
    
    for (i=1; i<nx; i++)
        R[i] = R[i-1]+ny;
    
    long Bt2 = nx*ny*ny*sizeof(double);
    
    R[0][0] = (double *)malloc(Bt2);
    
    if (!R[0][0])
        nrerror("allocation failure 3 in MATRIX3_RHO()");
    
    if (Flag_Initialize)
        memset(R[0][0], 0, Bt2);
    
    for (i=1; i<nx; i++) {
        R[i][0] = R[i-1][0]+ny*ny;
    }
    
    for (i=0; i<nx; i++) {
        for (j=1; j <ny; j++) {
            R[i][j] = R[i][j-1]+2*j;
        }
    }
    
    return R;
}


extern void FREE_MATRIX3_RHO_DB(double ***Rho){
    
    /******************************************************************************************
     Purpose:
     Free a complex rho matrix allocated by MATRIX3_RHO().
     Record of revisions:
     7th Jun. 2019
     Input parameters:
     ***Rho, the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/

    free(Rho[0][0]);
    
    free(Rho[0]);
    
    free(Rho);
    
    return;
}

extern double ***CUBE_DOUBLE(long Ni0, long Ni1, long Nj0, long Nj1, long Nk0, long Nk1, int Flag_Initialize){
    
    /******************************************************************************************
     Purpose:
     Allocate an complex double rho matrix m[i][j][k] with subscript range 0<=i<=n1,
     0<=j<=n2, -j<=k<=j.
     Record of revisions:
     7th Jun. 2019
     Input parameters:
     nh, the dimension of the rho matrix.
     Return:
     return the pointer.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/
    
    long Ni = Ni1-Ni0+1, Nj = Nj1-Nj0+1, Nk = Nk1-Nk0+1, i, j;
    double ***C;
    
    long Byte = Ni*sizeof(double **);
    C = (double ***)malloc(Byte);
    if (!C){
        nrerror("allocation failure 1 in CUBE_DOUBLE()");
    }
    C -= Ni0;
    
    Byte = (Ni*Nj)*sizeof(double * );
    C[Ni0] = (double **) malloc(Byte);
    if (!C[Ni0]){
        nrerror("allocation failure 2 in CUBE_DOUBLE()");
    }
    C[Ni0] -= Nj0;
    for(i=Ni0+1; i<=Ni1; i++){
        C[i] = C[i-1]+Nj;
    }
    
    Byte = (Ni*Nj*Nk)*sizeof(double);
    C[Ni0][Nj0] = (double *) malloc(Byte);
    if (!C[Ni0][Nj0]){
        nrerror("allocation failure 3 in CUBE_DOUBLE()");
    }
    C[Ni0][Nj0] -= Nk0;
    
    if (Flag_Initialize)
        memset(C[Ni0][Nj0], 0, Byte);
    
    for(i=Ni0+1; i<=Ni1; i++){
        C[i][Nj0] = C[i-1][Nj0]+Nj*Nk;
    }
    
    for(i=Ni0; i<=Ni1; i++){
        for (j=Nj0+1; j<=Nj1; j++) {
            C[i][j] = C[i][j-1]+Nk;
        }
    }
    
    
    return C;
}

extern void FREE_CUBE_DOUBLE(double ***C, long Ni0, long Nj0, long Nk0){
    
    /******************************************************************************************
     Purpose:
     Free a char matrix allocated with MATRIX_CHAR().
     Record of revisions:
     20 Nov 2019
     Input parameters:
     **m, the pointer.
     nrl, the first row index of subscript range.
     ncl, the first column index of subscript range.
     Reference:
     Numerical recipes in C 2ed.
     ******************************************************************************************/
    
    free(C[Ni0][Nj0]+Nk0);
    free(C[Ni0]+Nj0);
    free((C+Ni0));
    
    return;
}

/**********************************************************************************************/
