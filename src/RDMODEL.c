
#include "RDMODEL.h"

/**********************************************************************************************/

extern void Compute_Ion(STRUCT_ATOM Atom, STRUCT_PARA *Para){
    
    /******************************************************************************************
     Purpose:
     Compute the ion fraction.
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     Atom, a structure saved the atomic information.
     Output parameters:
     *Para, a structure saved the infromation of physical paramter.
     ******************************************************************************************/

    if (((*Para).Tempture<Atom.Ioniz[0][0]) || ((*Para).Tempture>Atom.Ioniz[0][Atom.Num_Ioniz-1])) {
        (*Para).Ion = 0;
        return;
    }
    
    int i, indx = -1;
    double DT1 = 0, DT2 = 0;
    for (i=0; i<Atom.Num_Ioniz-1; i++) {
        if (((*Para).Tempture>=Atom.Ioniz[0][i]) && ((*Para).Tempture<=Atom.Ioniz[0][i+1])) {
            indx = i;
            DT1 = ((*Para).Tempture-Atom.Ioniz[0][i])/(Atom.Ioniz[0][i+1]-Atom.Ioniz[0][i]);
            DT2 = 1-DT1;
            break;
        }
    }
    (*Para).Ion = (DT2*Atom.Ioniz[1][i]+DT1*Atom.Ioniz[1][i+1])*(*Para).Hydrogen*Atom.Abundance;
    
    return;
}


extern void Compute_Para(STRUCT_CORONAL_MODEL Model, STRUCT_SPH_COORD Position, STRUCT_PARA *Para){
    
    /******************************************************************************************
     Purpose:
     Compute physical parameters at a position of the model.
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     Model, a structure saved the coronal model.
     Position, the position.
     Output parameters:
     *Para, a structure saved the infromation of physical paramter.
     ******************************************************************************************/

    if ((Position.R<Model.R[0]) || (Position.R >Model.R[Model.R_Len-1])) {
        fprintf(stderr, "R is out of the model range \n");
	    fprintf(stderr, "%e %e %e \n", Position.R, Model.R[0], Model.R[Model.R_Len-1]);
        exit(1);
    }
    if ((Position.Theta<Model.Theta[0]) || (Position.Theta>Model.Theta[Model.Theta_Len-1])) {
        fprintf(stderr, "Theta is out of the model range \n");
        fprintf(stderr, "%e %e %e \n",Position.Theta, Model.Theta[0], Model.Theta[Model.Theta_Len-1]);
        if ((Position.Theta<Model.Theta[0])) {
            Position.Theta = Model.Theta[0];
        }else if (Position.Theta>Model.Theta[Model.Theta_Len-1]){
            Position.Theta = Model.Theta[Model.Theta_Len-1];
        }
        //exit(1);
    }
    if ((Position.Phi<Model.Phi[0]) || (Position.Phi >Model.Phi[Model.Phi_Len-1])) {
        fprintf(stderr, "Phi is out of the model range \n");
        fprintf(stderr, "%e %e %e \n",Position.Phi, Model.Phi[0], Model.Phi[Model.Phi_Len-1]);
        if ((Position.Phi<Model.Phi[0])) {
            Position.Phi = Model.Phi[0];
        }else if (Position.Phi >Model.Phi[Model.Phi_Len-1]){
            Position.Phi  = Model.Phi[Model.Phi_Len-1];
        }
        //exit(1);
    }
    
    int i, indx_r = -1, indx_theta = -1, indx_phi = -1;
    double Scale_T= 1e6;
    
    for (i=0; i<Model.R_Len-1; i++) {
        if ((Position.R>=Model.R[i]) && (Position.R <=Model.R[i+1])) {
            indx_r = i;
            break;
        }
    }
    
    for (i=0; i<Model.Theta_Len-1; i++) {
        if ((Position.Theta>=Model.Theta[i]) && (Position.Theta <=Model.Theta[i+1])) {
            indx_theta = i;
            break;
        }
    }
    
    for (i=0; i<Model.Phi_Len-1; i++) {
        if (Position.Phi >= Model.Phi[i] && Position.Phi <= Model.Phi[i+1]) {
            indx_phi = i;
            break;
        }
    }
    
    double Ratio_R = (Position.R-Model.R[indx_r])/(Model.R[indx_r+1]-Model.R[indx_r]);
    double Ratio_Theta = (Position.Theta-Model.Theta[indx_theta])/(Model.Theta[indx_theta+1] \
                                                                   -Model.Theta[indx_theta]);
    double Ratio_Phi = (Position.Phi-Model.Phi[indx_phi])/(Model.Phi[indx_phi+1] \
                                                           -Model.Phi[indx_phi]);
    double Density = Interpolation_Linear_3D(Model.Density, indx_r, indx_theta, indx_phi, \
                                             Ratio_R, Ratio_Theta, Ratio_Phi);
    
    double H_Density = C_H_Ratio/C_mH;
    double He_Density = C_He_Ratio/C_mHe*2;
    (*Para).Hydrogen = Density;
    (*Para).Electron = (*Para).Hydrogen*(H_Density+He_Density)/H_Density;
    
    double T = Interpolation_Linear_3D(Model.T, indx_r, indx_theta, indx_phi, Ratio_R, \
                                       Ratio_Theta, Ratio_Phi);
    (*Para).Tempture = T*Scale_T;

    double Br = Interpolation_Linear_3D(Model.Br, indx_r, indx_theta, indx_phi, Ratio_R, \
                                        Ratio_Theta, Ratio_Phi);
    double Btheta = Interpolation_Linear_3D(Model.Btheta, indx_r, indx_theta, indx_phi, \
                                            Ratio_R, Ratio_Theta, Ratio_Phi);
    double BPhi = Interpolation_Linear_3D(Model.Bphi, indx_r, indx_theta, indx_phi, Ratio_R, \
                                          Ratio_Theta, Ratio_Phi);
    double B = sqrt(Br*Br+Btheta*Btheta+BPhi*BPhi);
    
    (*Para).Mag.B = B;
    (*Para).Mag.ThetaB = acos(Br/B);
    (*Para).Mag.PhiB = atan2(BPhi, Btheta);
    
    return;
}


extern void ReadModel(STRUCT_CORONAL_MODEL *Model, char Filename[]){
    
    /******************************************************************************************
     Purpose:
     Read the coronal model.
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     *Model_Path, a structure saved the paths of the files.
     Output parameters:
     *Model, a structure saved the coronal model.
     ******************************************************************************************/

    int indx_r, indx_theta, indx_phi;
    FILE *Fa, *Fb;
    char file[300];
    
    Fb = fopen(Filename, "r");
    fgets(file, 300, Fb);
    Trim(file, 3);
    
    Fa = fopen(file, "r");
    (*Model).R = (double *)malloc((*Model).R_Len*sizeof(double));
    for (indx_r=0; indx_r<(*Model).R_Len; indx_r++) {
        fscanf(Fa, "%le ",&(*Model).R[indx_r]);
    }
    fclose(Fa);
    
    fgets(file, 300, Fb);
    Trim(file, 3);
    Fa=fopen(file, "r");
    (*Model).Theta = (double *)malloc((*Model).Theta_Len*sizeof(double));
    for (indx_theta=0; indx_theta<(*Model).Theta_Len; indx_theta++) {
        fscanf(Fa, "%le ",&(*Model).Theta[indx_theta]);
    }
    fclose(Fa);
    
    fgets(file, 300, Fb);
    Trim(file, 3);
    Fa=fopen(file, "r");
    (*Model).Phi = (double *)malloc((*Model).Phi_Len*sizeof(double));
    for (indx_phi=0; indx_phi<(*Model).Phi_Len; indx_phi++) {
        fscanf(Fa, "%le ",&(*Model).Phi[indx_phi]);
    }
    fclose(Fa);
    
    fgets(file, 300, Fb);
    Trim(file, 3);
    Fa=fopen(file, "r");
    (*Model).Br = CUBE_DOUBLE(0, (*Model).R_Len-1, 0, (*Model).Theta_Len-1, 0, (*Model).Phi_Len, 0);
    for (indx_phi=0; indx_phi<(*Model).Phi_Len; indx_phi++) {
        for (indx_theta=0; indx_theta<(*Model).Theta_Len; indx_theta++) {
            for (indx_r=0; indx_r<(*Model).R_Len; indx_r++) {
                fscanf(Fa, "%le ",&(*Model).Br[indx_r][indx_theta][indx_phi]);
            }
        }
    }
    fclose(Fa);
    
    fgets(file, 300, Fb);
    Trim(file, 3);
    Fa=fopen(file, "r");
    (*Model).Btheta = CUBE_DOUBLE(0, (*Model).R_Len-1, 0, (*Model).Theta_Len-1, 0, (*Model).Phi_Len, 0);
    for (indx_phi=0; indx_phi<(*Model).Phi_Len; indx_phi++) {
        for (indx_theta=0; indx_theta<(*Model).Theta_Len; indx_theta++) {
            for (indx_r=0; indx_r<(*Model).R_Len; indx_r++) {
                fscanf(Fa, "%le ",&(*Model).Btheta[indx_r][indx_theta][indx_phi]);
            }
        }
    }
    fclose(Fa);
    
    fgets(file, 300, Fb);
    Trim(file, 3);
    Fa=fopen(file, "r");
    (*Model).Bphi = CUBE_DOUBLE(0, (*Model).R_Len-1, 0, (*Model).Theta_Len-1, 0, (*Model).Phi_Len, 0);
    for (indx_phi=0; indx_phi<(*Model).Phi_Len; indx_phi++) {
        for (indx_theta=0; indx_theta<(*Model).Theta_Len; indx_theta++) {
            for (indx_r=0; indx_r<(*Model).R_Len; indx_r++) {
                fscanf(Fa, "%le ",&(*Model).Bphi[indx_r][indx_theta][indx_phi]);
            }
        }
    }
    fclose(Fa);
    
    fgets(file, 300, Fb);
    Trim(file, 3);
    Fa=fopen(file, "r");
    (*Model).T = CUBE_DOUBLE(0, (*Model).R_Len-1, 0, (*Model).Theta_Len-1, 0, (*Model).Phi_Len, 0);
    for (indx_phi=0; indx_phi<(*Model).Phi_Len; indx_phi++) {
        for (indx_theta=0; indx_theta<(*Model).Theta_Len; indx_theta++) {
            for (indx_r=0; indx_r<(*Model).R_Len; indx_r++) {
                fscanf(Fa, "%le ",&(*Model).T[indx_r][indx_theta][indx_phi]);
            }
        }
    }
    fclose(Fa);
    
    fgets(file, 300, Fb);
    Trim(file, 3);
    Fa=fopen(file, "r");
    (*Model).Density = CUBE_DOUBLE(0, (*Model).R_Len-1, 0, (*Model).Theta_Len-1, 0, (*Model).Phi_Len, 0);
    for (indx_phi=0; indx_phi<(*Model).Phi_Len; indx_phi++) {
        for (indx_theta=0; indx_theta<(*Model).Theta_Len; indx_theta++) {
            for (indx_r=0; indx_r<(*Model).R_Len; indx_r++) {
                fscanf(Fa, "%le ",&(*Model).Density[indx_r][indx_theta][indx_phi]);
            }
        }
    }
    fclose(Fa);
    fclose(Fb);
    
    return;
}


extern void FreeModel(STRUCT_CORONAL_MODEL *Model, STRUCT_PARA *Para){
    
    /******************************************************************************************
     Purpose:
     Free the memory of the coronal model.
     Record of revisions:
     30 Nov. 2019
     Input parameters:
     *Model, a structure saved the coronal model.
     ******************************************************************************************/

    free((*Model).R);
    free((*Model).Theta);
    free((*Model).Phi);
    FREE_CUBE_DOUBLE((*Model).Br, 0, 0, 0);
    FREE_CUBE_DOUBLE((*Model).Btheta, 0, 0, 0);
    FREE_CUBE_DOUBLE((*Model).Bphi, 0, 0, 0);
    FREE_CUBE_DOUBLE((*Model).Density, 0, 0, 0);
    FREE_CUBE_DOUBLE((*Model).T, 0, 0, 0);
    free(Model);

    free(Para);

    return;
}


extern double Interpolation_Linear_3D(double ***Data, int i, int j, int k, \
                                      double ir, double jr, double kr){
    
    /******************************************************************************************
     Purpose:
     3D linear Interpolation.
     Record of revisions:
     26 Sep. 2019
     Input parameters:
     Data[][][], 3D data cube.
     i, j, k, the grid point close the interpolation position.
     ir, jr, kr, the ratios.
     return:
     return the interpolated value.
     ******************************************************************************************/

    double tmp1, tmp2, tmp3, tmp4;
    double Ratio_XC = 1-ir;
    double Ratio_YC = 1-jr;

    tmp1 = Ratio_XC*Data[i][j][k]+ir*Data[i+1][j][k];
    tmp2 = Ratio_XC*Data[i][j+1][k]+ir*Data[i+1][j+1][k];
    tmp3 = Ratio_YC*tmp1+jr*tmp2;
    tmp1 = Ratio_XC*Data[i][j][k+1]+ir*Data[i+1][j][k+1];
    tmp2 = Ratio_XC*Data[i][j+1][k+1]+ir*Data[i+1][j+1][k+1];
    tmp4 = Ratio_YC*tmp1+jr*tmp2;
    
    return (1-kr)*tmp3+kr*tmp4;
}

/**********************************************************************************************/
