
#ifndef RDINPUT_h
#define RDINPUT_h

/******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "STRING2.h"
#include "LL04.h"

/******************************************************************************************/

typedef struct{
    int Verbose;
    int Single_Point;
    int Flag_Collsion;
    int Flag_Fullcalculation;
    int Flag_Reduced;
    int Flag_Symmetry;
    int Thomson_scattering;
    double X[2], Z[2], dy, Rmax;
    int Nmax;
    double X_pos, Z_pos, Y_pos;
    char Atom_Path[300];
    char Model_Path[300];
    char Output_Path[300];
}STRUCT_INPUT;

/******************************************************************************************/

extern void REINPUT(char Filename[], STRUCT_INPUT *Input, STRUCT_OUTPUT *Output);

/******************************************************************************************/

#endif /* RDINPUT_h */
