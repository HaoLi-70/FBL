
#ifndef RDMODEL_h
#define RDMODEL_h

/******************************************************************************************/

#include <stdio.h>
#include "CONSTANT.h"
#include "LL04.h"
#include "STRING2.h"

/******************************************************************************************/

typedef struct{
    int R_Len, Theta_Len, Phi_Len;
    double *R, *Theta, *Phi;
    double ***T;
    double ***Density;
    double ***Br;
    double ***Btheta;
    double ***Bphi;
}STRUCT_CORONAL_MODEL;

typedef struct Struct_Coronal_Path{
    char *R_path;
    char *Theta_path;
    char *Phi_path;
    char *T_path;
    char *Density_path;
    char *Br_path;
    char *Btheta_path;
    char *Bphi_path;
}STRUCT_MODEL_PATH;

/******************************************************************************************/

extern void Compute_Ion(STRUCT_ATOM Atom, STRUCT_PARA *Para);

extern void Compute_Para(STRUCT_CORONAL_MODEL Model, STRUCT_SPH_COORD Position, STRUCT_PARA *Para);

extern void ReadModel(STRUCT_CORONAL_MODEL *Model, char Filename[]);

extern void FreeModel(STRUCT_CORONAL_MODEL *Model, STRUCT_PARA *Para);

extern double Interpolation_Linear_3D(double ***Data, int i, int j, int k, \
                                      double ir, double jr, double kr);

/******************************************************************************************/

#endif /* RDMODEL_h */
