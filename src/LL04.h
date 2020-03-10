
#ifndef LL04_h
#define LL04_h

/**********************************************************************************************/

#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "LU.h"
#include "CONSTANT.h"

/**********************************************************************************************/

typedef struct{
    int au, al;
    double w[3];
    int Flag_Mag_Diple;
    double Emcoefi[4];
    double ***Stokes;
    int Indx_Transition;
    double Kr, Kt;
}STRUCT_TRANSITION_OUT;

typedef struct{
    int Num_transition;
    STRUCT_TRANSITION_OUT *TR;
}STRUCT_OUTPUT;



typedef struct{
    double Energy;
    double J, L, S;
    double g;
    complex double **Rho;
}STRUCT_LEVEL;

typedef struct Struct_Transition{
    int au, al;
    double nu, lambda, Intensity;
    complex double **J_KQ, **J_KQ0;
    double u1, u2;
}STRUCT_TRANSITION;

typedef struct {
    char Element[10];
    double Mass, Abundance;
    int Num_Level, Num_Transition, Num_EQ;
    STRUCT_LEVEL *LV;
    double **A, **B;
    int Num_T;
    double *T;
    double ***Col_Strength;
    int Num_Ioniz;
    double **Ioniz;
    STRUCT_TRANSITION *TR;
}STRUCT_ATOM;


typedef struct{
    double B;
    double ThetaB;
    double PhiB;
}STRUCT_MAG;

typedef struct{
    double Tempture;
    double Ion;
    double Electron;
    double Velocity;
    double Hydrogen;
    STRUCT_MAG Mag;
}STRUCT_PARA;


typedef struct{
    double R;
    double Theta;
    double Phi;
}STRUCT_SPH_COORD;

/**********************************************************************************************/

extern double Geffect(double Gu, double Gl, double Ju, double Jl);

extern double Gfactor(double J, double L, double S);

extern double TJ(double J1, double J2, double J3, double M1, double M2, double M3);

extern double SJ(double J1, double J2, double J3, double J4, double J5, double J6);

extern double NJ(double J1, double J2, double J3, double J4, double J5, double J6, \
                 double J7, double J8, double J9);

extern double complex Djmn(double Alpha, double Beta, double Gamma, double J, \
                           double M, double N);

extern complex double TA(int a, double J, int K, int Q, int al, double Jl, int Kl, int Ql, \
                         STRUCT_ATOM atom);

extern double TE(int a, double J, int K, int Q, int au, double Ju, int Ku, int Qu, \
                 STRUCT_ATOM Atom);

extern complex double TS(int a, double J, int K, int Q, int au, float Ju, int Ku, int Qu, \
                         STRUCT_ATOM Atom);

extern complex double RA(int a, double J, int K, int Q, int Kp, int Qp, STRUCT_ATOM Atom);

extern double RE(int a, double J, int K, int Q, int Kp, int Qp, STRUCT_ATOM Atom);

extern complex double RS(int a, double J, int K, int Q, int Kp, int Qp, STRUCT_ATOM Atom);

extern void TP_tensor (complex double ***T);

extern void TQ_tensor (complex double ***T, double Theta, double Chi, double Gamma);

double Aniso(double u1, double u2, double r);

double Intensity(double u1, double u2, double r);

void Limb_Darkening(double Lambda, double *u1, double *u2);

extern void Incident_Tensor(STRUCT_ATOM Atom, STRUCT_SPH_COORD Position, \
                            STRUCT_MAG Mag, int Flag_symmetry);

void Rho_Compute(STRUCT_ATOM Atom, STRUCT_SPH_COORD Position, STRUCT_PARA Para, \
                 double **Col_Rates, int Col_Flag, int Flag_Reduced);
    
void Init_Output(STRUCT_OUTPUT Output, STRUCT_ATOM Atom);

extern void Free_Output(STRUCT_OUTPUT *Output);

extern void EMCoefi_Inteall(STRUCT_ATOM Atom, STRUCT_PARA Para, STRUCT_OUTPUT output, \
                            int Flag_Cal, int Flag_Reduced);

extern void EMCoefi_Inte(STRUCT_ATOM Atom, complex double ***T_KQ, \
                         STRUCT_TRANSITION_OUT *Transition, int Flag_Cal, int Flag_Reduced);

extern void Thom_Scat_van(double r, double *PB, double q);

extern void Thom_Scat(double r, double *PB, double u1, double u2);

/**********************************************************************************************/

#endif /* LL04_h */
