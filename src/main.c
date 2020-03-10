
#include <stdio.h>
#include <math.h>
#include "CONSTANT.h"
#include "RDATOM.h"
#include "ALLOCATION.h"
#include "LL04.h"
#include "TIME_PRINT.h"
#include "RDMODEL.h"
#include "RDINPUT.h"

/******************************************************************************************/

#define Path_Input "input.dat"


/******************************************************************************************/

int main(int argc, const char * argv[]) {
    
    Time_Print();


    STRUCT_INPUT *Input;
    Input = (STRUCT_INPUT *)malloc(sizeof(STRUCT_INPUT));
    STRUCT_OUTPUT *Output;
    Output = (STRUCT_OUTPUT *)malloc(sizeof(STRUCT_OUTPUT));
    STRUCT_ATOM *Atom;
    Atom = (STRUCT_ATOM *) malloc(sizeof(STRUCT_ATOM));

    char *Filename=Path_Input;
    REINPUT(Filename, Input, Output);

    READATOM((*Input).Atom_Path, Atom);

    INIT_ATOM(Atom, (*Input).Flag_Reduced);

    Init_Output(*Output, *Atom);

    STRUCT_CORONAL_MODEL *Mod_Corona;
    Mod_Corona = (STRUCT_CORONAL_MODEL *)malloc(sizeof(STRUCT_CORONAL_MODEL));
    (*Mod_Corona).R_Len=150;
    (*Mod_Corona).Theta_Len=100;
    (*Mod_Corona).Phi_Len =181;

    ReadModel(Mod_Corona, (*Input).Model_Path);
    
    double **Col_Rates;
    Col_Rates=MATRIX_DOUBLE(0, (*Atom).Num_Level-1, 0, (*Atom).Num_Level-1, 1);
    
    STRUCT_PARA *Para;
    Para=(STRUCT_PARA *)malloc(sizeof(STRUCT_PARA));
    STRUCT_SPH_COORD Position;
    int i;
    double PB[2], u1, u2;
    double X, Y, Z, R;
    double Kr, Kt, cont2, COnST;
    
    Time_Print();

    if ((*Input).Single_Point) {
        
        X=(*Input).X_pos;
        Z=(*Input).Z_pos;
        Y=(*Input).Y_pos;
        
        R=sqrt(X*X+Y*Y+Z*Z);
        
        Position.R=R;
        Position.Theta = acos(Z / R);
        Position.Phi=atan2(Y, X);
        
        fprintf(stderr, "%e %e \n", Position.Theta, Position.Phi);
        
        Compute_Para(*Mod_Corona, Position, Para);
        
        Compute_Ion(*Atom, Para);
        
        Collisional_Rates(Atom, (*Para).Tempture, Col_Rates);
        
        Incident_Tensor((*Atom), Position, (*Para).Mag, (*Input).Flag_Symmetry);
        
        Rho_Compute(*Atom, Position, *Para, Col_Rates, (*Input).Flag_Collsion, \
                    (*Input).Flag_Reduced);
        
        EMCoefi_Inteall(*Atom, *Para, *Output, (*Input).Flag_Fullcalculation, \
                        (*Input).Flag_Reduced);

        for (i=0; i<(*Output).Num_transition; i++) {
            fprintf(stderr, "\nTransition i: upper level = %d, lower level = %d \n", \
                    (*Output).TR[i].au, (*Output).TR[i].al);
            fprintf(stderr, "I = %e, Q = %e, U = %e \n",(*Output).TR[i].Emcoefi[0]*(*Para).Ion, \
                    (*Output).TR[i].Emcoefi[1]*(*Para).Ion, (*Output).TR[i].Emcoefi[2]*(*Para).Ion);
            fprintf(stderr, "Q/I = %.4f, U/I = %.4f \n",(*Output).TR[i].Emcoefi[1] \
                    /(*Output).TR[i].Emcoefi[0],(*Output).TR[i].Emcoefi[2] \
                    /(*Output).TR[i].Emcoefi[0]);
            if((*Input).Thomson_scattering){
                u1 = (*Atom).TR[(*Output).TR[i].Indx_Transition].u1;
                u2 = (*Atom).TR[(*Output).TR[i].Indx_Transition].u2;
                Thom_Scat(R, PB, u1, u2);

                cont2 = C_c/(*Atom).TR[(*Output).TR[i].Indx_Transition].lambda \
                /(*Atom).TR[(*Output).TR[i].Indx_Transition].lambda \
                *(*Atom).TR[(*Output).TR[i].Indx_Transition].Intensity*(*Para).Electron;
                COnST = 3./8.*C_sigmaT*1e-10;

                Kr = ((1-(X*X+Z*Z)/R/R)*PB[0]+(X*X+Z*Z)/R/R*PB[1])*COnST*cont2;
                Kt = PB[0]*COnST*cont2;
                
                fprintf(stderr,"Kr = %e , Kt = %e \n",Kr,Kt);
                fprintf(stderr,"Thomson scattering I = %e , Q = %e franction = %e \n",Kr+Kt,Kt-Kr,(Kt-Kr)/(Kr+Kt));
            }
        }
        
        
    }else{
        
        int ix, iz, iy, Nstep, j;
        double weight;
        COnST = 3./8.*C_Solarradus*C_sigmaT*1e-10;
        
        for (j=0; j<(*Output).Num_transition; j++) {
            (*Output).TR[j].Stokes=CUBE_DOUBLE(1, (*Input).Nmax, 1, (*Input).Nmax, 0, 4, 1);
        }
        
        for (ix = 1; ix <= (*Input).Nmax; ix++) {
            X = (*Input).X[0]+((*Input).X[1]-(*Input).X[0])*(ix-1.)/((*Input).Nmax-1.);
            for (iz = 1; iz <= (*Input).Nmax; iz++) {
                Z = (*Input).Z[0]+((*Input).Z[1]-(*Input).Z[0])*(iz-1.)/((*Input).Nmax-1.);
            
                if ((*Input).Verbose) {
                    fprintf(stderr, "ix = %d, iz = %d, X= %e, Z= %e \n", ix, iz, X, Z);
                }
                if (X*X+Z*Z>1.05) {
                    
                    Kr = 0;
                    Kt = 0;
                    Nstep = (int)(sqrt((*Input).Rmax*(*Input).Rmax-X*X-Z*Z)/(*Input).dy);
                    
                    for(iy = -Nstep; iy <= Nstep; iy++){
                        Y = iy*(*Input).dy;
                        R = sqrt(X*X+Y*Y+Z*Z);
                        if (iy == -Nstep||iy == Nstep){
                            weight = 0.5;
                        }else{
                            weight = 1.;
                        }
                        
                        Position.R=R;
                        Position.Theta = acos(Z / R);
                        Position.Phi=atan2(Y, X);
                        if (Position.Phi<0) {
                            Position.Phi = C_Pi*2+Position.Phi;
                        }
                        
                        Compute_Para(*Mod_Corona, Position, Para);
                        
                        Collisional_Rates(Atom, (*Para).Tempture, Col_Rates);
                        
                        Incident_Tensor((*Atom), Position, (*Para).Mag, \
                                        (*Input).Flag_Symmetry);
                        
                        Rho_Compute(*Atom, Position, *Para, Col_Rates, \
                                    (*Input).Flag_Collsion, (*Input).Flag_Reduced);
                        
                        Compute_Ion(*Atom, Para);
                        
                        EMCoefi_Inteall(*Atom, *Para, *Output, (*Input).Flag_Fullcalculation, \
                                        (*Input).Flag_Reduced);
                        
                        for (j=0; j<(*Output).Num_transition; j++) {
                            
                            for (i=0; i<3; i++){
                                (*Output).TR[j].Stokes[ix][iz][i] += (*Output).TR[j].Emcoefi[i] \
                                *(*Para).Ion*(*Input).dy*weight*C_Solarradus;
                            }
                            if((*Input).Thomson_scattering){
                                u1 = (*Atom).TR[(*Output).TR[j].Indx_Transition].u1;
                                u2 = (*Atom).TR[(*Output).TR[j].Indx_Transition].u2;
                                Thom_Scat(R, PB, u1, u2);
        //the factor c/lambda/lambda is due to dnv = c/lambda/lambda*dlambda
                                cont2 = C_c/(*Atom).TR[(*Output).TR[j].Indx_Transition].lambda \
                                    /(*Atom).TR[(*Output).TR[j].Indx_Transition].lambda \
                                    *(*Atom).TR[(*Output).TR[j].Indx_Transition].Intensity \
                                    *(*Para).Electron*weight*(*Input).dy;
                                
                                Kr = ((1-(X*X+Z*Z)/R/R)*PB[0]+(X*X+Z*Z)/R/R*PB[1]);
                                Kt = PB[0];
                                (*Output).TR[j].Stokes[ix][iz][3] += (Kr + Kt)*COnST*cont2;
                                
                                (*Output).TR[j].Stokes[ix][iz][4] += (Kr - Kt)*COnST*cont2;
                               
        
                            }
                        }
                    }
                
                }else{
                    for (j=0; j<(*Output).Num_transition; j++) {
                        (*Output).TR[j].Stokes[ix][iz][0] += 1e-30;
                    }
                }
            }
        }
        Time_Print();

        FILE *Fa;
        Fa = fopen((*Input).Output_Path, "w");
        for (ix = 1; ix <= (*Input).Nmax; ix++) {
            for (iz = (*Input).Nmax; iz >= 1; iz--) {
            //for (iz = 1; iz <= (*Input).Nmax; iz++) {
                for (j=0; j<(*Output).Num_transition; j++) {
                    for (i=0; i<3; i++){
                        fprintf(Fa, "%e ",(*Output).TR[j].Stokes[ix][iz][i]);
                    }
                    if((*Input).Thomson_scattering){
                        for (i=3; i<5; i++){
                            fprintf(Fa, "%e ",(*Output).TR[j].Stokes[ix][iz][i]);
                        }
                    }
                }
                fprintf(Fa, "\n");
            }
        }
        for (i=0; i<(*Output).Num_transition; i++) {
            FREE_CUBE_DOUBLE((*Output).TR[i].Stokes, 1, 1, 0);
        }
        fclose(Fa);
    }
    
    Free_Output(Output);
    free(Input);
    FreeModel(Mod_Corona, Para);
    FREEATOM(Atom, Col_Rates);
    Time_Print();
    
    return 0;
}

/******************************************************************************************/

