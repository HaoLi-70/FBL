
#ifndef RDATOM_h
#define RDATOM_h

/******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "STRING2.h"
#include "ALLOCATION.h"
#include "LL04.h"

/******************************************************************************************/


extern void Collisional_Rates(STRUCT_ATOM *atom, double T, double **Col_Rates);

extern void INIT_ATOM(STRUCT_ATOM *atom, int Flag);

extern void READATOM(char filename[], STRUCT_ATOM *atom);

extern void FREEATOM(STRUCT_ATOM *atom, double **Col_Rates);

/******************************************************************************************/

#endif /* RDATOM_h */
