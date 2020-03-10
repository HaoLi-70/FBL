
#include "RDATOM.h"

/**********************************************************************************************/

static double Planck(double nu, double T);

/**********************************************************************************************/

static double Planck(double nu, double T){
    
    /******************************************************************************************
     Purpose:
     Compute the Planck function.
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     nu, the frequency.
     T, the temperature.
     return:
     return the intensity at frequency nu for temperature T.
     ******************************************************************************************/
    
    double Intensity = (2*C_h*nu*nu*nu/C_c/C_c)/(exp(C_h*nu/C_Kb/T)-1);
    return Intensity;
}


extern void Collisional_Rates(STRUCT_ATOM *atom, double T, double **Col_Rates){
    
    /******************************************************************************************
     Purpose:
     Interpolate a collisional rates matrix for a temperature T.
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     *Atom, a structure saved the atomic information.
     T, the temperature.
     Output parameters:
     Col_Rates[][], a matrix saved the collisional rates;
     ******************************************************************************************/
    
    double T_log = log10(T);
    
    int i, j, indx = 0;
    double DT1 = 0, DT2 = 0;
    double Strength_tmp, tmp;

    if ( (T_log < (*atom).T[0]) || (T_log > (*atom).T[(*atom).Num_T-1])) {
        for (i=0; i<(*atom).Num_Level; i++) {
            for (j=0; j<(*atom).Num_Level; j++) {
                Col_Rates[i][j] = 0;
            }
        }
    }else{
        indx = -1;
        for (i=0; i<(*atom).Num_T-2; i++) {
            if(T_log>=(*atom).T[i] && T_log<(*atom).T[i+1]){
                indx = i;
                DT1 = (T_log-(*atom).T[i])/((*atom).T[i+1]-(*atom).T[i]);
                DT2 = 1-DT1;
                break;
            }
        }
        
        for (i=0; i<(*atom).Num_Level; i++) {
            for (j=i+1; j<(*atom).Num_Level; j++) {
                Strength_tmp = (*atom).Col_Strength[i][j][indx]*DT2+(*atom).Col_Strength[i][j][indx+1]*DT1;
                tmp = 8.63e-6*Strength_tmp/sqrt(T);
                Col_Rates[j][i] = tmp/(2*(*atom).LV[j].J+1);
                Col_Rates[i][j] = tmp/(2*(*atom).LV[i].J+1)*exp(((*atom).LV[i].Energy-(*atom).LV[j].Energy)*C_h*C_c/T/C_Kb);
            }
        }
    }
    
    return;
}


extern void INIT_ATOM(STRUCT_ATOM *Atom, int Flag){
    
    /******************************************************************************************
     Purpose:
     Initialize the atom information (B matrix, limb darkening coefficients etc.).
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     *Atom, a structure saved the atomic information.
     Output parameters:
     *Atom, a structure saved the atomic information.
     ******************************************************************************************/
    
    (*Atom).TR = (STRUCT_TRANSITION *)malloc((*Atom).Num_Transition*sizeof(STRUCT_TRANSITION));
    (*Atom).B = MATRIX_DOUBLE(0, (*Atom).Num_Level-1, 0, (*Atom).Num_Level-1, 1);

    int i, j, indx = 0;
    (*Atom).Num_EQ = 0;
    double T = 5800;
    for (i=0; i<(*Atom).Num_Level; i++) {
        
        if(Flag){
            (*Atom).Num_EQ += (int)((*Atom).LV[i].J+1);
        }else{
            (*Atom).Num_EQ += (int)(((*Atom).LV[i].J*2+1)*((*Atom).LV[i].J*2+1));
        }
        (*Atom).LV[i].Rho = MATRIX_RHO((int)(2*(*Atom).LV[i].J), 1);
        (*Atom).LV[i].g = Gfactor((*Atom).LV[i].J, (*Atom).LV[i].L, (*Atom).LV[i].S);
        
        for (j=i+1; j<(*Atom).Num_Level; j++) {
            if((*Atom).A[j][i]>0){
                (*Atom).TR[indx].au = j;
                (*Atom).TR[indx].al = i;
                (*Atom).TR[indx].lambda = 1/((*Atom).LV[j].Energy-(*Atom).LV[i].Energy);
                (*Atom).TR[indx].nu = C_c/(*Atom).TR[indx].lambda;
                (*Atom).TR[indx].Intensity = Planck((*Atom).TR[indx].nu, T);
                (*Atom).TR[indx].J_KQ = MATRIX_RHO(2, 1);
                (*Atom).TR[indx].J_KQ0 = MATRIX_RHO(2, 1);
                (*Atom).B[j][i] = 1/(2*C_h*C_c)*(*Atom).A[j][i]*(*Atom).TR[indx].lambda \
                    *(*Atom).TR[indx].lambda*(*Atom).TR[indx].lambda;
                (*Atom).B[i][j] = (*Atom).B[j][i]*(2*(*Atom).LV[j].J+1)/(2*(*Atom).LV[i].J+1);
                Limb_Darkening((*Atom).TR[indx].lambda, &(*Atom).TR[indx].u1, &(*Atom).TR[indx].u2);
                indx++;
            }
        }
    }
    return;
}


extern void FREEATOM(STRUCT_ATOM *Atom, double **Col_Rates){
    
    /******************************************************************************************
     Purpose:
     Free the memory of atomic structure and collisional rates.
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     *Atom, a structure saved the atomic information.
     Col_Rates, the pointer of the matrix.
     ******************************************************************************************/
    
    int i;
    for (i=0; i<(*Atom).Num_Transition; i++) {
        FREE_MATRIX_RHO((*Atom).TR[i].J_KQ);
        FREE_MATRIX_RHO((*Atom).TR[i].J_KQ0);
    }
    free((*Atom).TR);
    
    for (i=0; i<(*Atom).Num_Level; i++) {
        FREE_MATRIX_RHO((*Atom).LV[i].Rho);
    }
    free((*Atom).LV);
    
    FREE_VECTOR_DOUBLE((*Atom).T, 0);
    
    FREE_MATRIX_DOUBLE((*Atom).B, 0, 0);
    
    FREE_MATRIX_DOUBLE((*Atom).A, 0, 0);
    
    FREE_MATRIX_DOUBLE((*Atom).Ioniz, 0, 0);
    
    FREE_CUBE_DOUBLE((*Atom).Col_Strength, 0, 0, 0);
    
    FREE_MATRIX_DOUBLE(Col_Rates, 0, 0);
    
    free(Atom);
    
    return;
}


extern void READATOM(char filename[], STRUCT_ATOM *Atom){
    
    /******************************************************************************************
     Purpose:
     Read the atomic information from a file used for the forbidden line calculation.
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     Filename[], the input file.
     Output parameters:
     *Atom, a structure saved the atomic information.
     ******************************************************************************************/
    
    FILE *fa;
    fa = fopen(filename, "r");
    
    char lines[500], lines_tmp[500], *p, parameter[30];
    long len_tot;
    int flag = 0, indx = 0, levelindx[2], len, i, j;
    double T_tmp[200], Ratio[200], tmp[5], tmp2 = 0;

    while (fgets(lines, 300, fa) != NULL) {
        
        Trim(lines, 3);
        len_tot = strlen(lines);
        
        if(len_tot>0){
            
            if( (lines[0] != '#') && (lines[0] != '!')){
                
                switch (flag) {
                    case 0:
                        memmove((*Atom).Element,lines,10);
                        flag++;
                        break;
                        
                    case 1:
                        (*Atom).Mass = atof(lines);
                        flag++;
                        break;
                        
                    case 2:
                        (*Atom).Abundance = atof(lines);
                        flag++;
                        break;
                        
                    case 3:
                        (*Atom).Num_Level = atoi(lines);
                        flag++;
                        break;
                        
                    case 4:
                        indx = 0;
                        (*Atom).LV = (STRUCT_LEVEL *) malloc((*Atom).Num_Level*sizeof(STRUCT_LEVEL));
                        do {
                            if(len_tot>0){
                                p = lines;
                                for (i=0; i<5; i++) {
                                    memmove(lines_tmp, p, len_tot);
                                    lines_tmp[len_tot]='\0';
                                    len = Indx_Char(lines_tmp, ',', 1);
                                    if (len >1 ) {
                                        memmove(parameter, p, len-1);
                                        parameter[len-1] = '\0';
                                        p += len;
                                        len_tot -= len;
                                        tmp[i] = atof(parameter);
                                    }else{
                                        tmp[i] = atof(p);
                                    }
                                }
                                
                                (*Atom).LV[indx].Energy = tmp[1]*100;
                                (*Atom).LV[indx].J = tmp[2];
                                (*Atom).LV[indx].L = tmp[3];
                                (*Atom).LV[indx].S = tmp[4];
                                indx++;
                            }
                            
                            if(indx < (*Atom).Num_Level){
                                fgets(lines, 300, fa);
                                Trim(lines, 3);
                                len_tot = strlen(lines);
                            }
                        } while (indx < (*Atom).Num_Level);
                        flag++;
                        break;
                        
                    case 5:
                        (*Atom).A = MATRIX_DOUBLE(0, (*Atom).Num_Level-1, 0, (*Atom).Num_Level-1, 1);
                        (*Atom).Num_Transition = 0;
                        do {
                            if(len_tot>0){
                                (*Atom).Num_Transition++;
                                p = lines;
                                for (i=0; i<3; i++) {
                                    memmove(lines_tmp, p, len_tot);
                                    lines_tmp[len_tot] = '\0';
                                    len = Indx_Char(lines_tmp, ',', 1);
                                    if (len > 1 ) {
                                        memmove(parameter, p, len-1);
                                        parameter[len-1] = '\0';
                                        p += len;
                                        len_tot -= len;
                                        memmove(lines_tmp, p, len_tot);
                                        lines_tmp[len_tot] = '\0';
                                        levelindx[i] = atoi(parameter);
                                    }else{
                                        tmp2 = atof(p);
                                    }
                                }
                                (*Atom).A[levelindx[0]][levelindx[1]] = tmp2;
                            }
                            if (fgets(lines, 300, fa) == NULL){
                                break;
                            }else if ((lines[0] == '#') || (lines[0] == '!')){
                                break;
                            }else{
                                Trim(lines, 3);
                                len_tot = strlen(lines);
                            }
                        } while (1);
                        flag++;
                        break;
                        
                    case 6:
                        p = lines;
                        
                        for (i=0; ; i++) {
                            memmove(lines_tmp, p, len_tot);
                            lines_tmp[len_tot] = '\0';
                            len = Indx_Char(lines_tmp, ',', 1);
                            if (len > 1 ) {
                                
                                memmove(parameter, p, len-1);
                                parameter[len-1] = '\0';
                                p += len;
                                len_tot -= len;
                                memmove(lines_tmp, p, len_tot);
                                
                                lines_tmp[len_tot] = '\0';
                                T_tmp[i] = atof(parameter);
                            }else{
                                T_tmp[i] = atof(p);
                                break;
                            }
                        }
                        (*Atom).Num_T = i+1;
                        (*Atom).T = VECTOR_DOUBLE(0, i, 0);
                        for (j=0; j<(*Atom).Num_T; j++) {
                            (*Atom).T[j] = T_tmp[j];
                        }
                        flag++;
                        break;
                        
                    case 7:
                        (*Atom).Col_Strength = CUBE_DOUBLE(0, (*Atom).Num_Level-1, 0, \
                            (*Atom).Num_Level-1, 0, (*Atom).Num_T-1, 1);
                        do {
                            if(len_tot>0){
                                p = lines;
                                for (i=0,j=-2; j<(*Atom).Num_T; i++,j++) {
                                    memmove(lines_tmp, p, len_tot);
                                    lines_tmp[len_tot] = '\0';
                                    len = Indx_Char(lines_tmp, ',', 1);
                                    if (len >1 ) {
                                        memmove(parameter, p, len-1);
                                        parameter[len-1] = '\0';
                                        p += len;
                                        len_tot -= len;
                                        memmove(lines_tmp, p, len_tot);
                                        
                                        lines_tmp[len_tot] = '\0';
                                        if (j<0) {
                                            levelindx[i] = atoi(parameter);
                                        }else{
                                            (*Atom).Col_Strength[levelindx[0]][levelindx[1]][j]=atof(parameter);
                                        }
                                    }else{
                                        (*Atom).Col_Strength[levelindx[0]][levelindx[1]][j]=atof(p);
                                    }
                                }
                            }
                            if (fgets(lines, 300, fa) == NULL){
                                break;
                            }else if ((lines[0] == '#') || (lines[0] == '!')){
                                break;
                            }else{
                                Trim(lines, 3);
                                len_tot = strlen(lines);
                            }
                        } while (1);
                        flag++;
                        break;
                        
                    case 8:
                        indx = 0;
                        do {
                            if(len_tot>0){
                                p = lines;
                                len = Indx_Char(lines, ',', 1);
                                memmove(parameter, p, len-1);
                                parameter[len-1] = '\0';
                                p += len;
                                
                                T_tmp[indx] = atof(parameter);
                                Ratio[indx] = atof(p);
                                indx++;
                            }
                            if (fgets(lines, 300, fa) == NULL){
                                break;
                            }else if ((lines[0] =='#') || (lines[0] == '!')){
                                break;
                            }else{
                                Trim(lines, 3);
                                len_tot = strlen(lines);
                            }
                        } while (1);
                        (*Atom).Num_Ioniz = indx;
                        (*Atom).Ioniz = MATRIX_DOUBLE(0, 1, 0, indx-1, 0);
                        for (i=0; i<indx; i++) {
                            (*Atom).Ioniz[0][i] = T_tmp[i];
                            (*Atom).Ioniz[1][i] = Ratio[i];
                        }
                        flag++;
                        break;
                        
                    default:
                        break;
                }
            }
        }
    }
    fclose(fa);
    
    return;
}

/**********************************************************************************************/
