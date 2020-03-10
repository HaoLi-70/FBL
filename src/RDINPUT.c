
#include "RDINPUT.h"

/**********************************************************************************************/

extern void REINPUT(char Filename[], STRUCT_INPUT *Input, STRUCT_OUTPUT *Output){
    
    /******************************************************************************************
     Purpose:
     Read the input file for the forbidden line calculation.
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     Filename[], the input file.
     Output parameters:
     *Input, a structure saved the input information.
     *Output, a structure saved the output information.
     ******************************************************************************************/
    
    FILE *fa;
    fa = fopen(Filename, "r");
    
    char lines[300], lines_tmp[300], key[50], parameter[50], *p;
    int len, i, indx, level[50];
    long len_tot;
    
    (*Input).Verbose = 1;
    
    while (fgets(lines, 300, fa) != NULL) {
        Trim(lines, 3);
        len_tot = strlen(lines);
        
        if(len_tot > 0){
            len = Indx_Char(lines, '=', 1);
            len_tot -= len;
            
            if( (lines[0] != '#') && (lines[0] != '!') ){
                String_Copy(key, lines, len-1, 1);
                
                if ( strcmp(key, "verbose") == 0){
                    p = lines+len;
                    String_Copy(parameter, p, len_tot, 1);
                    
                    if (strcmp(parameter,"YES") == 0) {
                        (*Input).Verbose = 1;
                    }else{
                        (*Input).Verbose = 0;
                    }
                    if ((*Input).Verbose) {
                        fprintf(stderr, "\n Verbose : YES \n");
                    }else{
                        fprintf(stderr, "\n Verbose : NO \n");
                    }
                }else if ( strcmp(key, "single_point") == 0){
                    p = lines+len;
                    String_Copy(parameter, p, len_tot, 1);
                    
                    if (strcmp(parameter,"YES")==0) {
                        (*Input).Single_Point=1;
                    }else{
                        (*Input).Single_Point=0;
                    }
                    if ((*Input).Verbose) {
                        if ((*Input).Single_Point) {
                            fprintf(stderr, "\n Single Point Calculation: YES \n");
                        }else{
                            fprintf(stderr, "\n Single Point Calculation: No \n");
                        }
                    }
                }else if ( strcmp(key, "collision") == 0){
                    p = lines+len;
                    String_Copy(parameter, p, len_tot, 1);
                    
                    if (strcmp(parameter,"YES") == 0) {
                        (*Input).Flag_Collsion=1;
                    }else{
                        (*Input).Flag_Collsion = 0;
                    }
                    if ((*Input).Verbose) {
                        if ((*Input).Flag_Collsion) {
                            fprintf(stderr, "\n Collision Calculation: YES \n");
                            fprintf(stderr, "\n Warning: The multipolar components of collisional rates \n");
                            fprintf(stderr, "are calculated accoriding to the equations in LL04 appendx A1.4 \n");
                        }else{
                            fprintf(stderr, "\n Collision Calculation: No \n");
                        }
                    }
                }else if ( strcmp(key, "full_calculation") == 0){
                    p = lines+len;
                    String_Copy(parameter, p, len_tot, 1);
                    if (strcmp(parameter,"YES") == 0) {
                        (*Input).Flag_Fullcalculation = 1;
                    }else{
                        (*Input).Flag_Fullcalculation = 0;
                    }
                    if ((*Input).Verbose) {
                        if ((*Input).Flag_Fullcalculation) {
                            fprintf(stderr, "\n Full Calculation of the Emission Coefficients: YES \n");
                        }else{
                            fprintf(stderr, "\n Full Calculation of the Emission Coefficients: No \n");
                        }
                    }
                }else if ( strcmp(key, "symmetry ") == 0){
                    p = lines+len;
                    String_Copy(parameter, p, len_tot, 1);
                    if (strcmp(parameter,"YES") == 0) {
                        (*Input).Flag_Symmetry = 1;
                    }else{
                        (*Input).Flag_Symmetry = 0;
                    }
                    if ((*Input).Verbose) {
                        if ((*Input).Flag_Symmetry) {
                            fprintf(stderr, "\n Cylindrical symmetry : YES \n");
                        }else{
                            fprintf(stderr, "\n Cylindrical symmetry : No \n");
                        }
                    }
                }else if ( strcmp(key, "reduced") == 0){
                    p = lines+len;
                    String_Copy(parameter, p, len_tot, 1);
                    if (strcmp(parameter,"YES") == 0) {
                        (*Input).Flag_Reduced = 1;
                    }else{
                        (*Input).Flag_Reduced = 0;
                    }
                    if ((*Input).Verbose) {
                        if ((*Input).Flag_Reduced) {
                            fprintf(stderr, "\n Redue the SEEs : YES \n");
                        }else{
                            fprintf(stderr, "\n Redue the SEEs : No \n");
                        }
                    }
                }else if ( strcmp(key, "Thomson_scattering") == 0){
                    p = lines+len;
                    String_Copy(parameter, p, len_tot, 1);
                    if (strcmp(parameter,"YES") == 0) {
                        (*Input).Thomson_scattering = 1;
                    }else{
                        (*Input).Thomson_scattering = 0;
                    }
                    if ((*Input).Verbose) {
                        if ((*Input).Thomson_scattering) {
                            fprintf(stderr, "\n Compute Thomson scattering : YES \n");
                        }else{
                            fprintf(stderr, "\n Compute Thomson scattering : No \n");
                        }
                    }
                }else if (strcmp(key, "atom_path") == 0){
                    p = lines+len;
                    String_Copy((*Input).Atom_Path, p, len_tot, 1);
                    
                    if ((*Input).Verbose) {
                        fprintf(stderr, "\n Atom file : %s \n",(*Input).Atom_Path);
                    }
                }else if (strcmp(key, "output_path") == 0){
                    p = lines+len;
                    String_Copy((*Input).Output_Path, p, len_tot, 1);
                    
                    if ((*Input).Verbose) {
                        fprintf(stderr, "\n Output file : %s \n",(*Input).Output_Path);
                    }
                }else if (strcmp(key, "model_path") == 0){
                    p = lines+len;
                    String_Copy((*Input).Model_Path, p, len_tot, 1);
                    
                    if ((*Input).Verbose) {
                        fprintf(stderr, "\n Model file : %s \n",(*Input).Model_Path);
                    }
                }else if (strcmp(key, "x") == 0){
                    p = lines+len;
                    String_Copy(lines_tmp, p, len_tot, 0);
                    len = Indx_Char(lines_tmp, ',', 1);
                    String_Copy(parameter, p, len-1, 1);
                    p += len;
                    (*Input).X[0] = atof(parameter);
                    (*Input).X[1] = atof(p);
                    if ((*Input).Verbose&&!(*Input).Single_Point) {
                        fprintf(stderr, "\n X ranges from %e to %e \n",(*Input).X[0],(*Input).X[1]);
                    }
                }else if (strcmp(key, "z") == 0){
                    p = lines+len;
                    String_Copy(lines_tmp, p, len_tot, 0);
                    len = Indx_Char(lines_tmp, ',', 1);
                    String_Copy(parameter, p, len-1, 1);
                    p += len;
                    (*Input).Z[0] = atof(parameter);
                    (*Input).Z[1] = atof(p);
                    if ((*Input).Verbose&&!(*Input).Single_Point) {
                        fprintf(stderr, "\n Z ranges from %e to %e \n",(*Input).Z[0],(*Input).Z[1]);
                    }
                }else if (strcmp(key, "dy") == 0){
                    p = lines+len;
                    (*Input).dy = atof(p);
                    if ((*Input).Verbose&&!(*Input).Single_Point) {
                        fprintf(stderr, "\n dy : %e Rs \n",(*Input).dy);
                    }
                }else if (strcmp(key, "rmax") == 0){
                    p = lines+len;
                    (*Input).Rmax = atof(p);
                    if ((*Input).Verbose&&!(*Input).Single_Point) {
                        fprintf(stderr, "\n Max of the radius : %e Rs \n",(*Input).Rmax);
                    }
                }else if (strcmp(key, "nmax") == 0){
                    p = lines+len;
                    (*Input).Nmax = atoi(p);
                    if ((*Input).Verbose&&!(*Input).Single_Point) {
                        fprintf(stderr, "\n Nmax : %d \n",(*Input).Nmax);
                    }
                }else if (strcmp(key, "x_pos") == 0){
                    p = lines+len;
                    (*Input).X_pos = atof(p);
                    if ((*Input).Verbose&&(*Input).Single_Point) {
                        fprintf(stderr, "\n X Position : %e \n",(*Input).X_pos);
                    }
                }else if (strcmp(key, "z_pos") == 0){
                    p = lines+len;
                    (*Input).Z_pos = atof(p);
                    if ((*Input).Verbose&&(*Input).Single_Point) {
                        fprintf(stderr, "\n Z Position : %e \n",(*Input).Z_pos);
                    }
                }else if (strcmp(key, "y_pos") == 0){
                    p = lines+len;
                    (*Input).Y_pos = atof(p);
                    if ((*Input).Verbose&&(*Input).Single_Point) {
                        fprintf(stderr, "\n Y Position : %e \n",(*Input).Y_pos);
                    }
                }else if (strcmp(key, "upperlevel") == 0){
                    p = lines+len;
                    String_Copy(lines_tmp, p, len_tot, 0);
                    len = Indx_Char(lines_tmp, ',', 1);
                    indx = 0;
                    while (len > 1) {
                        String_Copy(parameter, p, len-1, 1);
                        p += len;
                        len_tot -= len;
                        String_Copy(lines_tmp, p, len_tot, 0);
                        level[indx] = atoi(parameter);
                        indx++;
                        len = Indx_Char(lines_tmp, ',', 1);
                    }
                    level[indx] = atoi(lines_tmp);
                    indx++;
                    (*Output).Num_transition = indx;
                    (*Output).TR = (STRUCT_TRANSITION_OUT *)malloc((*Output).Num_transition*sizeof(STRUCT_TRANSITION_OUT));
                    for (i=0; i<indx; i++) {
                        (*Output).TR[i].au = level[i];
                    }
                    if ((*Input).Verbose){
                        fprintf(stderr, "\n Output transition: \n");
                        fprintf(stderr, "Upper level index : ");
                        for (i=0; i<indx; i++) {
                            fprintf(stderr, " %d ",(*Output).TR[i].au);
                        }
                        fprintf(stderr, "\n");
                    }
                }else if (strcmp(key, "lowerlevel") == 0){
                    p = lines+len;
                    String_Copy(lines_tmp, p, len_tot, 0);
                    len = Indx_Char(lines_tmp, ',', 1);
                    indx = 0;
                    while (len > 1) {
                        String_Copy(parameter, p, len-1, 1);
                        p += len;
                        len_tot -= len;
                        String_Copy(lines_tmp, p, len_tot, 0);
                        level[indx] = atoi(parameter);
                        indx++;
                        len = Indx_Char(lines_tmp, ',', 1);
                    }
                    level[indx] = atoi(lines_tmp);
                    indx++;
                    if ((*Output).Num_transition != indx) {
                        fprintf(stderr, "\n Numbers of upper and lower levels do not match. ");
                        exit(1);
                    }
                    for (i=0; i<indx; i++) {
                        (*Output).TR[i].al = level[i];
                    }
                    if ((*Input).Verbose){
                        fprintf(stderr, "\n Output transition: \n");
                        fprintf(stderr, "Lower level index : ");
                        for (i=0; i<indx; i++) {
                            fprintf(stderr, " %d ",(*Output).TR[i].al);
                        }
                        fprintf(stderr, "\n");
                    }
                }else if (strcmp(key, "magneticdipole") == 0){
                    p = lines+len;
                    String_Copy(lines_tmp, p, len_tot, 0);
                    len = Indx_Char(lines_tmp, ',', 1);
                    indx = 0;
                    while (len > 1) {
                        String_Copy(parameter, p, len-1, 1);
                        if (strcmp(parameter,"YES") == 0) {
                            (*Output).TR[indx].Flag_Mag_Diple = 1;
                        }else{
                            (*Output).TR[indx].Flag_Mag_Diple = 0;
                        }
                        indx++;
                        p += len;
                        len_tot -= len;
                        String_Copy(lines_tmp, p, len_tot, 0);
                        len = Indx_Char(lines_tmp, ',', 1);
                    }
                    Trim(lines_tmp, 3);
                    if (strcmp(lines_tmp,"YES") == 0) {
                        (*Output).TR[indx].Flag_Mag_Diple = 1;
                    }else{
                        (*Output).TR[indx].Flag_Mag_Diple = 0;
                    }
                    
                    if ((*Input).Verbose){
                        fprintf(stderr, "\n Magnetic Dipole Transition: \n");
                        for (i=0; i<=indx; i++) {
                            if ((*Output).TR[i].Flag_Mag_Diple) {
                                fprintf(stderr, " Yes ");
                            }else{
                                fprintf(stderr, " No ");
                            }
                        }
                        fprintf(stderr, "\n");
                    }
                }else{
                    if ((*Input).Verbose) {
                        fprintf(stderr, "Warning : Neglect key words %s \n",key);
                        fprintf(stderr, "The whole line is %s \n",lines);
                    }
                }
                
            }
        }
    }
    fclose(fa);
    
    return;
}

/**********************************************************************************************/

