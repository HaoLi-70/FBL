
#include "LL04.h"

/**********************************************************************************************/

static double Factorial(int n);

static void SEE(STRUCT_ATOM Atom, complex double **SEE_Matrix, double **InCol_Rates, \
                int InCol_Flag, double Density_Electron, STRUCT_MAG Mag);

static void SEE2(STRUCT_ATOM Atom, complex double **SEE_Matrix, double **InCol_Rates, \
                 int InCol_Flag, double Density_Electron, STRUCT_MAG Mag);

/**********************************************************************************************/

static double Factorial(int n){
    
    /******************************************************************************************
     Purpose:
     Computes the factorial of an integer.
     Record of revisions:
     5 Jun 2019
     Input parameters:
     n, the integer.
     Return:
     return the factorial.
     ******************************************************************************************/
    
    if (n<=0) {
        return 1;
    } else {
        int i;
        double F = 1;
        for (i=n; i>0; i--) {
            F = F*i;
        }
        return F;
    }
}


extern double Geffect(double Gu, double Gl, double Ju, double Jl){
    
    /***********************************************************************************
     Purpose:
     Computes the effective Lande factor.
     Modified:
     25 Apr. 2018
     Input parameters:
     Gu, The lande factor of upper level.
     Gl, The lande factor of lower level.
     Ju, The total angular momentum of upper level.
     Jl, The total angular momentum of lower level.
     Output parameters:
     Geffect, The effect Lande factor.
     References:
     LL04 Chapter 3, Equation 3.44.
     ***********************************************************************************/
    
    double Geffect;
    
    Geffect = 0.5*(Gu+Gl)+0.25*(Gu-Gl)*(Ju*(Ju+1)-Jl*(Jl+1));
    
    return Geffect;
}


extern double Gfactor(double J, double L, double S){
    
    /***********************************************************************************
     Purpose:
     Computes the Lande factor.
     Modified:
     25 Apr. 2018
     Input parameters:
     J, The total angular momentum of the electronic cloud.
     L, The total orbital angular momentum of the electronic cloud.
     S, The total spin of the electronic cloud.
     Output parameters:
     Gfactor, The Lande factor.
     Method:
     L-S coupling asumption, and if the number is 0, return 0 for convernience.
     ***********************************************************************************/
    
    if (J==0)
        return 0;
    
    double g;
    g = 1.0+(J*(J+1.0)-L*(L+1.0)+S*(S+1.0))/(2.0*J*(J+1.0));
    
    return g;
}


extern double TJ(double J1, double J2, double J3, double M1, double M2, double M3){
    
    /***********************************************************************************
     Purpose:
     Computes the 3j symbol.
     Modified:
     25 Apr. 2018
     Input parameters:
     J1, J2, J3, The total angular momentum of the electronic cloud.
     M1, M2, M3, The magnetic quantum number.
     Output parameters:
     F, The 3j symbol.
     References:
     LL04 Chapter 2, Page 36, Equation 2.19 and Page 38, Equation 2.22.
     ***********************************************************************************/
    
    if((M1+M2+M3)!=0)
        return 0;
    
    if((J1+J2)<J3||fabs(J1-J2)>J3)
        return 0;
    
    if(fabs(M1)>J1)
        return 0;
    
    if(fabs(M2)>J2)
        return 0;
    
    if(fabs(M3)>J3)
        return 0;
    
    if(fmod((J1+J2+J3)*2.0,2.0)!=0.0)
        return 0;
    
    if(fmod((J1-M1)*2.0,2.0)!=0.0)
        return 0;
    
    if(fmod((J2-M2)*2.0,2.0)!=0.0)
        return 0;
    
    double a, c, d, e = 0, F, t, t1 = 0, t2;
    
    a = pow(-1,J1-J2-M3);
    c = sqrt(1.0*Factorial(J1+J2-J3)*Factorial(J1-J2+J3)*Factorial(-J1+J2+J3) \
             /Factorial(J1+J2+J3+1));
    d = sqrt(1.0*Factorial(J1+M1)*Factorial(J1-M1)*Factorial(J2+M2)*Factorial(J2-M2) \
             *Factorial(J3+M3)*Factorial(J3-M3));
    
    t1 = 0;
    if((J2-J3-M1)>t1)
        t1 = J2-J3-M1;
    
    if((J1+M2-J3)>t1)
        t1 = J1+M2-J3;
    
    t2 = J1+J2-J3;
    
    if((J1-M1)<t2)
        t2 = J1-M1;
    
    if ((J2+M2)<t2)
        t2 = J2+M2;
    
    for (t=t1; t<=t2; t++)
        e = e+pow(-1,t)/(1.0*Factorial(t)*Factorial(J1+J2-J3-t)*Factorial(J1-M1-t) \
                         *Factorial(J2+M2-t)*Factorial(J3-J2+M1+t)*Factorial(J3-J1-M2+t));
    
    F = a*c*d*e;
    return F;
}


extern double SJ(double J1, double J2, double J3, double J4, double J5, double J6){
    
    /***********************************************************************************
     Purpose:
     Computes the 6j symbol.
     Modified:
     25 Apr. 2018
     Input parameters:
     J1, J2, J3, J4, J5, J6, The total angular momentum of the electronic cloud.
     Output parameters:
     C, The 6j symbol.
     References:
     LL04 Chapter 2, Page 42, Equation 2.35.
     ***********************************************************************************/
    
    if ((J1+J2)<J3||fabs(J1-J2)>J3)
        return 0;
    
    if ((J1+J5)<J6||fabs(J1-J5)>J6)
        return 0;
    
    if ((J4+J2)<J6||fabs(J4-J2)>J6)
        return 0;
    
    if ((J4+J5)<J3||fabs(J4-J5)>J3)
        return 0;
    
    if(fmodf((J1+J2+J3)*2.0,2.0)!=0)
        return 0;
    
    if(fmodf((J1+J5+J6)*2.0,2.0)!=0)
        return 0;
    
    if(fmodf((J4+J2+J6)*2.0,2.0)!=0)
        return 0;
    
    if(fmodf((J4+J5+J3)*2.0,2.0)!=0)
        return 0;
    
    double a, b = 0, C, t1, t2, t;
    
    a = sqrt(1.0*Factorial(J1+J2-J3)*Factorial(J1-J2+J3)*Factorial(-J1+J2+J3) \
            /Factorial(J1+J2+J3+1.0))*sqrt(1.0*Factorial(J1+J5-J6)*Factorial(J1-J5+J6) \
            *Factorial(-J1+J5+J6)/Factorial(J1+J5+J6+1.0))*sqrt(1.0*Factorial(J4+J2-J6) \
            *Factorial(J4-J2+J6)*Factorial(-J4+J2+J6)/Factorial(J4+J2+J6+1.0)) \
            *sqrt(1.0*Factorial(J4+J5-J3)*Factorial(J4-J5+J3)*Factorial(-J4+J5+J3) \
            /Factorial(J4+J5+J3+1.0));
    
    t1 = J1+J2+J3;
    
    if((J1+J5+J6)>t1)
        t1 = J1+J5+J6;
    
    if((J4+J2+J6)>t1)
        t1 = J4+J2+J6;
    
    if((J4+J5+J3)>t1)
        t1=J4+J5+J3;
    
    t = t1;
    t2 = J1+J2+J4+J5;
    
    if((J2+J3+J5+J6)<t2)
        t2 = J2+J3+J5+J6;
    
    if((J1+J3+J4+J6)<t2)
        t2 = J1+J3+J4+J6;
    
    for (t=t1; t<=t2; t++)
        b = b+pow(-1, t)*Factorial(t+1)/Factorial(t-J1-J2-J3)/Factorial(t-J1-J5-J6) \
            /Factorial(t-J4-J2-J6)/Factorial(t-J4-J5-J3)/(Factorial(J1+J2+J4+J5-t) \
            *Factorial(J2+J3+J5+J6-t)*Factorial(J1+J3+J4+J6-t));
    
    C = a*b;
    return C;
}


extern double NJ(double J1, double J2, double J3, double J4, double J5, double J6, \
                 double J7, double J8, double J9){
    
    /***********************************************************************************
     Purpose:
     Computes the 9j symbol.
     Modified:
     25 Apr. 2018
     Input parameters:
     J1, J2, J3, J4, J5, J6, J7, J8, J9, The total angular momentum of the electronic cloud.
     Output parameters:
     C, The 9j symbol.
     References:
     LL04 Chapter 2, Page 47, Equation 2.48.
     ***********************************************************************************/
    
    if ((J1+J2)<J3||fabs(J1-J2)>J3)
        return 0;
    
    if ((J4+J5)<J6||fabs(J4-J5)>J6)
        return 0;
    
    if ((J7+J8)<J9||fabs(J7-J8)>J9)
        return 0;
    
    if ((J1+J4)<J7||fabs(J1-J4)>J7)
        return 0;
    
    if ((J2+J5)<J8||fabs(J2-J5)>J8)
        return 0;
    
    if ((J3+J6)<J9||fabs(J3-J6)>J9)
        return 0;
    
    if(fmod((J1+J2+J3)*2.0,2.0)!=0)
        return 0;
    
    if(fmod((J4+J5+J6)*2.0,2.0)!=0)
        return 0;
    
    if(fmod((J7+J8+J9)*2.0,2.0)!=0)
        return 0;
    
    if(fmod((J1+J4+J7)*2.0,2.0)!=0)
        return 0;
    
    if(fmod((J2+J5+J8)*2.0,2.0)!=0)
        return 0;
    
    if(fmod((J3+J6+J9)*2.0,2.0)!=0)
        return 0;
    
    double t1, t2, t, a=0.0;
    t2 = J1+J9;
    
    if((t2>J2+J6))
        t2 = J2+J6;
    
    if(t2>(J4+J8))
        t2 = J4+J8;
    
    t1 = fabs(J1-J9);
    
    if(t1<fabs(J2-J6))
        t1 = fabs(J2-J6);
    
    if(t1<fabs(J8-J4))
        t1 = fabs(J8-J4);
    
    t = t1;
    
    for (t=t1; t<=t2; t++)
        a = a+pow(-1, 2*t)*(2*t+1)*SJ(J1, J9, t, J8, J4, J7) \
            *SJ(J2, J6, t, J4, J8, J5)*SJ(J1, J9, t, J6, J2, J3);
    
    return a;
}


extern double complex Djmn(double Alpha, double Beta, double Gamma, double J, double M, \
                           double N){
    
    /***********************************************************************************
     Purpose:
     Compute the element of the reduced rotation matrices.
     Modified:
     1 Dec. 2019
     Input parameters:
     Alpha, Beta, Gamma, the angles.
     J, M, N, quantum numbers.
     Return:
     the element of the reduced rotation matrices .
     References:
     LL04 Chapter 2, page 54 equation 2.68 and page 56 equation 2.69, rotation through
     an angle Alpha, Beta, and Gamma about the Z-axis, Y-axis, and Z-axis.
     ***********************************************************************************/

    double a, b = 0;
    double complex c, Djmn;
    int t, t1, t2;
    a = sqrt(Factorial((int)(J+M))*Factorial((int)(J-M))*Factorial((int)(J+N)) \
           *Factorial((int)(J-N)));
    
    t1 = (int)(M-N)<=0?0:(int)(M-N);
    
    t2 = (int)(J+M)<=(int)(J-N)?(int)(J+M):(int)(J-N);
    
    for (t=t1; t<=t2; t++) {
        b = b+pow(-1,t)*pow(cos(Beta/2),(2*J+M-N-2*t))*pow(sin(Beta/2),(2*t-M+N)) \
            /Factorial((int)(J+M-t))/Factorial((int)(J-N-t))/Factorial(t)/Factorial((int)(t+N-M));
    }
    
    c = cos(M*Alpha+N*Gamma)-sin(M*Alpha+N*Gamma)*I;
    Djmn = a*b*c;
    
    return Djmn;
}


extern complex double TA(int a, double J, int K, int Q, int al, double Jl, int Kl, int Ql, \
                         STRUCT_ATOM Atom){
    
    /***********************************************************************************
     Purpose:
     Transfer rates (TA in 7.14a).
     Modified:
     1 Dec. 2019
     Input parameters:
     a, J, K, Q, the quantum numbers of current level.
     al, Jl, Kl, Ql, the quantum numbers of lower level.
     Atom, a structure saved atomic information.
     Return:
     Transfer rates (TA).
     References:
     LL04 Chapter 7, page 287 equation 7.14a.
     ***********************************************************************************/

    if (Atom.B[al][a]<=0)
        return 0;

    int Kr, Qr, indx = -1, i;
    complex double ta=0;

    for (i=0; i<Atom.Num_Transition; i++) {
        if (a==Atom.TR[i].au && al==Atom.TR[i].al) {
            indx = i;
            break;
        }
    }
    
    if (indx>=0 && indx <Atom.Num_Transition) {
        for(Kr=0;Kr<=2;Kr++){
            for(Qr=-Kr;Qr<=Kr;Qr++){
                ta += (2*Jl+1)*Atom.B[al][a]*sqrt(3.0*(2*K+1)*(2*Kl+1)*(2*Kr+1)) \
                    *pow(-1,Kl+Ql)*NJ(J,Jl,1,J,Jl,1,K,Kl,Kr)*TJ(K,Kl,Kr,-Q,Ql,-Qr) \
                    *Atom.TR[indx].J_KQ[Kr][Qr];
            }
        }
        ta *= Atom.TR[indx].Intensity;
    }
    return ta;
}


extern double TE(int a, double J, int K, int Q, int au, double Ju, int Ku, int Qu, \
                 STRUCT_ATOM Atom){
    
    /***********************************************************************************
     Purpose:
     Transfer rates (TE in 7.14b).
     Modified:
     1 Dec. 2019
     Input parameters:
     a, J, K, Q, the quantum numbers of current level.
     au, Ju, Ku, Qu, the quantum numbers of upper level.
     Atom, a structure saved atomic information.
     Return:
     Transfer rates (TE).
     References:
     LL04 Chapter 7, page 287 equation 7.14b.
     ***********************************************************************************/
    
    if (Atom.A[au][a]<=0) return 0;
    double te;
    
    te = (K==Ku)*(Q==Qu)*(2*Ju+1)*Atom.A[au][a]*pow(-1,1+(int)(J+Ju)+K)*SJ(Ju,Ju,K,J,J,1);
    return te;
}


extern complex double TS(int a, double J, int K, int Q, int au, float Ju, int Ku, int Qu, \
                         STRUCT_ATOM Atom){
    
    /***********************************************************************************
     Purpose:
     Transfer rates (TS in 7.14c).
     Modified:
     1 Dec. 2019
     Input parameters:
     a, J, K, Q, the quantum numbers of current level.
     au, Ju, Ku, Qu, the quantum numbers of upper level.
     Atom, a structure saved atomic information.
     Return:
     Transfer rates (TS).
     References:
     LL04 Chapter 7, page 287 equation 7.14c.
     ***********************************************************************************/
    
    if (Atom.B[au][a]<=0)
        return 0;
    
    int Kr, Qr, i, indx = -1;
    double ts = 0;

    for (i=0; i<Atom.Num_Transition; i++) {
        if (au==Atom.TR[i].au && a==Atom.TR[i].al) {
            indx = i;
            break;
        }
    }
    
    if (indx>=0 && indx <Atom.Num_Transition) {
        for(Kr=0;Kr<=2;Kr++){
            for(Qr=-Kr;Qr<=Kr;Qr++){
                ts += (2*Ju+1)*Atom.B[au][a]*sqrt(3*(2*K+1)*(2*Ku+1)*(2*Kr+1)) \
                    *pow(-1,Kr+Ku+Qu)*NJ(J,Ju,1,J,Ju,1,K,Ku,Kr)*TJ(K,Ku,Kr,-Q,Qu,-Qr) \
                    *Atom.TR[indx].J_KQ[Kr][Qr];
            }
        }
        ts *= Atom.TR[indx].Intensity;
    }
    return ts;
}


extern complex double RA(int a, double J, int K, int Q, int Kp, int Qp, STRUCT_ATOM Atom){
    
    /***********************************************************************************
     Purpose:
     Transfer rates (RA in 7.14d).
     Modified:
     1 Dec. 2019
     Input parameters:
     a, J, K, Q, Kp, Qp, the quantum numbers of current level.
     Atom, a structure saved atomic information.
     Return:
     Transfer rates (RA).
     References:
     LL04 Chapter 7, page 288 equation 7.14d.
     ***********************************************************************************/
    
    if (a>=Atom.Num_Level-1||a<0)
        return 0;
    
    int Kr, Qr, au, indx, i;
    complex double ra = 0, temp = 0;
    double Ju;

    for (au=a+1; au<Atom.Num_Level; au++) {
        if (Atom.B[a][au]>0) {
            indx = -1;
            for (i=0; i<Atom.Num_Transition; i++) {
                if (au==Atom.TR[i].au && a==Atom.TR[i].al) {
                    indx = i;
                    break;
                }
            }
            if (indx>=0 && indx<Atom.Num_Transition) {
                Ju = Atom.LV[au].J;
                temp = 0;
                for(Kr=0.0;Kr<=2;Kr++){
                    for(Qr=-Kr;Qr<=Kr;Qr++){
                        temp += (2*J+1)*Atom.B[a][au]*sqrt(3*(2*K+1)*(2*Kp+1)*(2*Kr+1)) \
                            *pow(-1,(int)(1+(int)(Ju-J)+Kr+Qp))*SJ(K,Kp,Kr,J,J,J) \
                            *SJ(1,1,Kr,J,J,Ju)*TJ(K,Kp,Kr,Q,-Qp,Qr)*(1+pow(-1,K+Kp+Kr)) \
                            /2*Atom.TR[indx].J_KQ[Kr][Qr];
                    }
                }
                ra += temp*Atom.TR[indx].Intensity;
            }
        }
    }
    return ra;
}


extern double RE(int a, double J, int K, int Q, int Kp, int Qp, STRUCT_ATOM Atom){
    
    /***********************************************************************************
     Purpose:
     Transfer rates (RE in 7.14e).
     Modified:
     1 Dec. 2019
     Input parameters:
     a, J, K, Q, Kp, Qp, the quantum numbers of current level.
     Atom, a structure saved atomic information.
     Return:
     Transfer rates (RE).
     References:
     LL04 Chapter 7, page 288 equation 7.14e.
     ***********************************************************************************/
    
    if (a<=0)
        return 0;
    
    double re = 0;
    int al;
    
    if (Kp==K&&Qp==Q){
        for (al=0; al<a; al++)
            re += Atom.A[a][al];
    }else{
      return 0;
    }
    return re;
}

extern complex double RS(int a, double J, int K, int Q, int Kp, int Qp, STRUCT_ATOM Atom){
    
    /***********************************************************************************
     Purpose:
     Transfer rates (RS in 7.14f).
     Modified:
     1 Dec. 2019
     Input parameters:
     a, J, K, Q, Kp, Qp, the quantum numbers of current level.
     Atom, a structure saved atomic information.
     Return:
     Transfer rates (RS).
     References:
     LL04 Chapter 7, page 288 equation 7.14f.
     ***********************************************************************************/
    
    if (a>=Atom.Num_Level-1||a<0)
        return 0;
    
    int Kr, Qr, al, indx, i;
    complex double rs = 0, temp = 0;
    double Jl;
    
    for (al = 0; al < a; al++) {
        if (Atom.A[a][al] > 0) {
            indx=-1;
            for (i=0; i<Atom.Num_Transition; i++) {
                if (a==Atom.TR[i].au && al==Atom.TR[i].al) {
                    indx = i;
                    break;
                }
            }
            if (indx>=0 && indx <Atom.Num_Transition) {
                Jl = Atom.LV[al].J;
                temp = 0;
                for(Kr = 0; Kr <= 2; Kr += 2){
                    for(Qr = -Kr; Qr <= Kr; Qr++){
                        temp += (2*J+1)*Atom.B[a][al]*sqrt(3.*(2.*K+1.)*(2.*Kp+1.) \
                            *(2.*Kr+1.))*pow(-1,1+(int)(Jl-J)+Qp)*SJ(K, Kp, Kr, J, J, J) \
                            *SJ(1, 1, Kr, J, J, Jl)*TJ(K, Kp, Kr, Q, -Qp, Qr) \
                            *(1+pow(-1,K+Kp+Kr))/2.*Atom.TR[indx].J_KQ[Kr][Qr];
                    }
                }
                rs += temp*Atom.TR[indx].Intensity;
            }
        }
    }
    return rs;
}


extern void TP_tensor (complex double ***T){
    
    /***********************************************************************************
     Purpose:
     Geometry tensor TKP.
     Modified:
     1 Dec. 2019
     Output parameters:
     T[i][k][p], the tensors.
     References:
     LL04 Chapter 5, page 210 table 5.5.
     ***********************************************************************************/
    
    int i, K, P;
    for (i = 0; i < 4; i++) {
        for (K = 0; K < 3; K++) {
            for (P = -K; P <= K; P++) {
                T[i][K][P] = 0;
            }
        }
    }
    T[0][0][0] = 1;
    T[0][2][0] = 1/sqrt(2.);
    T[1][2][-2] = -sqrt(3.)/2;
    T[1][2][2] = -sqrt(3.)/2;
    T[2][2][-2] = sqrt(3.)/2*I;
    T[2][2][2] = -sqrt(3.)/2*I;
    T[3][1][0] = sqrt(3.)/2;
    
    return;
}


extern void TQ_tensor (complex double ***T, double Theta, double Chi, double Gamma){
    
    /***********************************************************************************
     Purpose:
     Geometry tensor TKQ.
     Modified:
     1 Dec. 2019
     Input parameters:
     Theta, Chi, Gamma, the angles in radians.
     Output parameters:
     T[i][k][p], the tensors.
     References:
     LL04 Chapter 5, page 210 table 5.6.
     ***********************************************************************************/
    
    // i=0
    T[0][0][0] = 1.;
    T[0][1][-1] = 0;
    T[0][1][0] = 0;
    T[0][1][1] = 0;
    T[0][2][-2] = sqrt(3.)/4*sin(Theta)*sin(Theta)*cos(2*Chi)+sqrt(3.)/4*sin(Theta) \
        *sin(Theta)*sin(2*Chi)*I;
    T[0][2][-1] = sqrt(3.)/2*sin(Theta)*cos(Theta)*cos(Chi)-sqrt(3.)/2*sin(Theta) \
        *cos(Theta)*sin(Chi)*I;
    T[0][2][0] = 1/2.0/sqrt(2.)*(3*cos(Theta)*cos(Theta)-1);
    T[0][2][1] = -sqrt(3.)/2*sin(Theta)*cos(Theta)*cos(Chi)-sqrt(3.)/2*sin(Theta) \
        *cos(Theta)*sin(Chi)*I;
    T[0][2][2] = sqrt(3.)/4*sin(Theta)*sin(Theta)*cos(2*Chi)-sqrt(3.)/4*sin(Theta) \
        *sin(Theta)*sin(2*Chi)*I;
    
    // i=1
    T[1][0][0] = 0;
    T[1][1][-1] = 0;
    T[1][1][0] = 0;
    T[1][1][1] = 0;
    T[1][2][-2] = -sqrt(3.)/4*(cos(2*Gamma)*(1+cos(Theta)*cos(Theta))*cos(2*Chi)- \
        2*sin(2*Gamma)*cos(Theta)*sin(2*Chi))+sqrt(3.)/4*(cos(2*Gamma)*(1+cos(Theta) \
        *cos(Theta))*sin(2*Chi)+2*sin(2*Gamma)*cos(Theta)*cos(2*Chi))*I;
    T[1][2][-1] = sqrt(3.)/2*sin(Theta)*(cos(2*Gamma)*cos(Theta)*cos(Chi)-sin(2*Gamma) \
        *sin(Chi))+-sqrt(3.)/2*sin(Theta)*(cos(2*Gamma)*cos(Theta)*sin(Chi)+sin(2*Gamma) \
        *cos(Chi))*I;
    T[1][2][0] = -3/2.0/sqrt(2)*cos(2*Gamma)*sin(Theta)*sin(Theta);
    T[1][2][1] = -sqrt(3.)/2*sin(Theta)*(cos(2*Gamma)*cos(Theta)*cos(Chi)-sin(2*Gamma) \
        *sin(Chi))-sqrt(3.)/2*sin(Theta)*(cos(2*Gamma)*cos(Theta)*sin(Chi)+sin(2*Gamma) \
        *cos(Chi))*I;
    T[1][2][2] = -sqrt(3.)/4*(cos(2*Gamma)*(1+cos(Theta)*cos(Theta))*cos(2*Chi)- \
        2*sin(2*Gamma)*cos(Theta)*sin(2*Chi))-sqrt(3.)/4*(cos(2*Gamma)*(1+cos(Theta) \
        *cos(Theta))*sin(2*Chi)+2*sin(2*Gamma)*cos(Theta)*cos(2*Chi))*I;
    
    // i=2
    T[2][0][0] = 0;
    T[2][1][-1] = 0;
    T[2][1][0] = 0;
    T[2][1][1] = 0;
    T[2][2][-2] = sqrt(3.)/4*(sin(2*Gamma)*(1+cos(Theta)*cos(Theta))*cos(2*Chi)+ \
        2*cos(2*Gamma)*cos(Theta)*sin(2*Chi))-sqrt(3.)/4*(sin(2*Gamma)*(1+cos(Theta) \
        *cos(Theta))*sin(2*Chi)-2*cos(2*Gamma)*cos(Theta)*cos(2*Chi))*I;
    T[2][2][-1] = -sqrt(3.)/2*sin(Theta)*(sin(2*Gamma)*cos(Theta)*cos(Chi)+cos(2*Gamma) \
        *sin(Chi))+sqrt(3.)/2*sin(Theta)*(sin(2*Gamma)*cos(Theta)*sin(Chi)-cos(2*Gamma) \
        *cos(Chi))*I;
    T[2][2][0] = -3/2.0/sqrt(2.)*sin(2*Gamma)*sin(Theta)*sin(Theta);
    T[2][2][1] = sqrt(3.)/2*sin(Theta)*(sin(2*Gamma)*cos(Theta)*cos(Chi)+cos(2*Gamma) \
        *sin(Chi))+sqrt(3.)/2*sin(Theta)*(sin(2*Gamma)*cos(Theta)*sin(Chi)-cos(2*Gamma) \
        *cos(Chi))*I;
    T[2][2][2] = sqrt(3.)/4*(sin(2*Gamma)*(1+cos(Theta)*cos(Theta))*cos(2*Chi)+ \
        2*cos(2*Gamma)*cos(Theta)*sin(2*Chi))+sqrt(3.)/4*(sin(2*Gamma)*(1+cos(Theta) \
        *cos(Theta))*sin(2*Chi)-2*cos(2*Gamma)*cos(Theta)*cos(2*Chi))*I;
    
    // i=3
    T[3][0][0] = 0;
    T[3][1][-1] = sqrt(3.)/2*sin(Theta)*cos(Chi)-sqrt(3.)/2*sin(Theta)*sin(Chi)*I;
    T[3][1][0] = sqrt(3.)/sqrt(2.)*cos(Theta);
    T[3][1][1] = -sqrt(3.)/2*sin(Theta)*cos(Chi)-sqrt(3.)/2*sin(Theta)*sin(Chi)*I;
    T[3][2][-2] = 0;
    T[3][2][-1] = 0;
    T[3][2][0] = 0;
    T[3][2][1] = 0;
    T[3][2][2] = 0;
    
    return;
}


double Intensity(double u1, double u2, double r){
    
    /***********************************************************************************
     Purpose:
     Compute incident radiation tensor J00.
     Modified:
     1 Dec. 2019
     Input parameters:
     u1, u2, the limb darkening coefficients
     r, the radius.
     Output parameters:
     the radiation tensor J00.
     References:
     LL04 Chapter 12, page 675 Eqs. 12.34-12.37.
     ***********************************************************************************/
    
    double Sr, Cr, a[3];
    Sr = 1/r;
    Cr = sqrt(1-Sr*Sr);
    if (u1 ==0 && u2 ==0) {
        return 1/2.*(1-Cr);
    }
    
    a[0] = 1-Cr;
    a[1] = Cr-0.5-0.5*Cr*Cr/Sr*log((1+Sr)/Cr);
    a[2] = (Cr+2)*(Cr-1)/3/(Cr+1);
    
    double Jv=0.5*(a[0]+a[1]*u1+a[2]*u2);
    return Jv;
}


double Aniso(double u1, double u2, double r){
    
    /***********************************************************************************
     Purpose:
     Compute incident radiation tensor J20.
     Modified:
     1 Dec. 2019
     Input parameters:
     u1, u2, the limb darkening coefficients
     r, the radius.
     Output parameters:
     the radiation tensor J20.
     References:
     LL04 Chapter 12, page 675 Eqs. 12.34-12.37.
     ***********************************************************************************/
    
    double Sr, Cr, c[3];
    Sr=1/r;
    Cr=sqrt(1-Sr*Sr);
    if (u1 ==0 && u2 ==0) {
        return 1/4./sqrt(2.)*Cr*Sr*Sr;
    }
    c[0]=Cr*Sr*Sr;
    c[1]=1/8.*(8*Cr*Cr*Cr-3*Cr*Cr-8*Cr+2+(4-3*Cr*Cr)*Cr*Cr/Sr*log((1+Sr)/Cr));
    c[2]=(Cr-1)/15./(Cr+1)*(9*Cr*Cr*Cr+18*Cr*Cr+7*Cr-4);
    
    double Kn = (c[0]+c[1]*u1+c[2]*u2);//3Kv-Jv;
    return Kn/sqrt(2.)/4.;
}


void Limb_Darkening(double Lambda, double *u1, double *u2){
    
    /***********************************************************************************
     Purpose:
     Interpolate the limb darkening coefficients at a wavlegth of Lambda.
     Modified:
     1 Dec. 2019
     Input parameters:
     Lambda, the wavelength (m^-1).
     Output parameters:
     *u1, *u2, the limb darkening coefficients.
     References:
     Astrophysical Quantities 3rd ed. - C. Allen (Athlone Press, 1973)
     ***********************************************************************************/
    
    double Lambda_C[22] = {0.20,0.22,0.245,0.265,0.28,0.30,0.32,0.35,0.37,0.38,0.40,\
        0.45,0.50,0.55,0.60,0.80,1.0,1.5,2.0,3.0,5.0,10.0};
    double u1_C[22] = {0.12,-1.3,-0.1,-0.1,0.38,0.74,0.88,0.98,1.03,0.92,0.91,0.99,\
        0.97,0.92,0.88,0.73,0.64,0.57,0.48,0.35,0.22,0.15};
    double u2_C[22] = {0.33,1.6,0.85,0.90,0.57,0.20,0.03,-0.1,-0.16,-0.05,-0.05,\
        -0.17,-0.22,-0.23,-0.23,-0.22,-0.20,-0.21,-0.18,-0.12,-0.07,-0.07};
    
    int i, indx = -1;
    double lambda_tmp = Lambda*1e6;
    
    if (lambda_tmp<Lambda_C[0]||lambda_tmp>Lambda_C[21]) {
        *u1 = 0;
        *u2 = 0;
        return;
    }
    
    for (i=0; i<21; i++) {
        if (lambda_tmp>=Lambda_C[i]&&lambda_tmp<Lambda_C[i+1]) {
            indx = i;
            break;
        }
    }
    double l1 = (lambda_tmp-Lambda_C[indx])/(Lambda_C[indx+1]-Lambda_C[indx]);
    double l2 = (Lambda_C[indx+1]-lambda_tmp)/(Lambda_C[indx+1]-Lambda_C[indx]);
    
    *u1 = l1*u1_C[indx+1]+l2*u1_C[indx];
    *u2 = l1*u2_C[indx+1]+l2*u2_C[indx];
    return;
}


extern void Incident_Tensor(STRUCT_ATOM Atom, STRUCT_SPH_COORD Position, \
                            STRUCT_MAG Mag, int Flag_symmetry){
    
    /***********************************************************************************
     Purpose:
     Compute incident radiation in the reference of magnetic field for each transition.
     Modified:
     1 Dec. 2019
     Input parameters:
     Atom, a structure saved the atom information.
     Position, the position.
     Mag, a structure saved the magnetic field.
     Flag_symmetry, if Flag_symmetry, only J00 and J20 are not zero in the reference
     of local vertical.
     Output parameters:
     Atom, a structure saved the atom information.
     References:
     LL04.
     ***********************************************************************************/
    
    complex double **J_KQ0;
    J_KQ0 = MATRIX_RHO(2, 1);
    J_KQ0[0][0] = Intensity(0, 0, Position.R);
    J_KQ0[2][0] = Aniso(0, 0, Position.R);
    
    int i, j, k, m;
    double tmp;
    
    for (i=0; i<Atom.Num_Transition; i++) {
        if (Atom.TR[i].u1==0 && Atom.TR[i].u2==0) {
            Atom.TR[i].J_KQ0[0][0] = J_KQ0[0][0];
            Atom.TR[i].J_KQ0[2][0] = J_KQ0[2][0];
        }else{
            Atom.TR[i].J_KQ0[0][0] = Intensity(Atom.TR[i].u1, Atom.TR[i].u2, Position.R);
            Atom.TR[i].J_KQ0[2][0] = Aniso(Atom.TR[i].u1, Atom.TR[i].u2, Position.R);
        }
    }
    
    for (i=0; i<Atom.Num_Transition; i++) {
        for (j=-2; j<=2; j++) {
            Atom.TR[i].J_KQ[2][j] = 0;
        }
        Atom.TR[i].J_KQ[0][0] = Atom.TR[i].J_KQ0[0][0];
    }
    
    k=0;
    
    if(Flag_symmetry){
        for (i=-2; i<=2; i++) {
            for (j=-2; j<=2; j++) {
                tmp =Djmn(0, -Position.Theta, -Position.Phi, 2, 0, j) \
                *Djmn(Mag.PhiB, Mag.ThetaB, 0, 2, j, i);
                for (m=0; m<Atom.Num_Transition; m++) {
                    Atom.TR[m].J_KQ[2][i] +=  Atom.TR[m].J_KQ0[2][0]*tmp;
                }
            }
        }
    }else{
        for (i=-2; i<=2; i++) {
            for (j=-2; j<=2; j++) {
                for (k=-2; k<=2; k++) {
                    tmp =Djmn(0, -Position.Theta, -Position.Phi, 2, k, j) \
                        *Djmn(Mag.PhiB, Mag.ThetaB, 0, 2, j, i);
                    for (m=0; m<Atom.Num_Transition; m++) {
                        Atom.TR[m].J_KQ[2][i] +=  Atom.TR[m].J_KQ0[2][k]*tmp;
                    }
                }
            }
        }
    }
    FREE_MATRIX_RHO(J_KQ0);
    return;
}


static void SEE(STRUCT_ATOM Atom, complex double **SEE_Matrix, double **InCol_Rates, \
                        int InCol_Flag, double Density_Electron, STRUCT_MAG Mag){
    
    /***********************************************************************************
     Purpose:
     Compute statistical equilibrium equations coefficients for multi-level atom.
     Modified:
     1 Dec. 2019
     Input parameters:
     Atom, a structure saved the atom information.
     InCol_Rates, inelastic collisional rates.
     InCol_Flag, inelastic collision flag.
     Density_Electron, electron density.
     Mag, a structure saved the magnetic field.
     Output parameters:
     SEE_Matrix, the coefficient matrix.
     References:
     LL04, Chapter 7, Page 285, Eq 711.
     ***********************************************************************************/
    
    int i, j, K1, K2, Q1, Q2, indx, a = 0, b = 0;
    double tmp1 = 0, tmp2 = 0;
    
    for (i=0; i<Atom.Num_Level; i++) {
        for (K1=0; K1<=Atom.LV[i].J*2; K1++) {
            for (Q1=-K1; Q1<=K1; Q1++) {
                b = 0;
                for (j=0; j<Atom.Num_Level; j++) {
                    if (InCol_Flag > 0) {
                        tmp1 = SJ(Atom.LV[i].J, Atom.LV[i].J, 0, Atom.LV[j].J, \
                                  Atom.LV[j].J, 1);
                        
                        if (tmp1 != 0 ) {
                            tmp2 = pow(-1, K1)*SJ(Atom.LV[i].J, Atom.LV[i].J, K1, \
                                                  Atom.LV[j].J, Atom.LV[j].J, 1)/
                            SJ(Atom.LV[i].J, Atom.LV[i].J, 0, Atom.LV[j].J, \
                               Atom.LV[j].J, 1);
                           
                        }else{
                            if (K1==0) {
                                tmp2 = 1;
                            }else{
                                tmp2 = 1;
                            }
                        }
                       
                    }
                    for (K2=0; K2<=Atom.LV[j].J*2; K2++) {
                        for (Q2=-K2; Q2<=K2; Q2++) {
                            if (i>j) {
                                if (Atom.B[i][j]>0) {
                                    SEE_Matrix[a][b] = TA(i, Atom.LV[i].J, K1, Q1, \
                                        j, Atom.LV[j].J, K2, Q2, Atom);
                                }else{
                                    SEE_Matrix[a][b]=0;
                                }

                                if (InCol_Flag > 0 && K1==K2 && Q1==Q2) {
                                    if (tmp2 != 0) {
                                        SEE_Matrix[a][b] += sqrt((2*Atom.LV[j].J+1) \
    /(2*Atom.LV[i].J+1))*tmp2*InCol_Rates[j][i]*Density_Electron;
                                    }
                                }
                            }else if(i==j){
                                SEE_Matrix[a][b] = -RA(i, Atom.LV[i].J, K1, Q1, K2, Q2, Atom) \
    -RE(i, Atom.LV[i].J, K1, Q1, K2, Q2, Atom)-RE(i, Atom.LV[i].J, K1, Q1, K2, Q2, Atom) \
    -RS(i, Atom.LV[i].J, K1, Q1, K2, Q2, Atom);
                                if (a==b && Q1!=0) {
                                    SEE_Matrix[a][b] -= 2.*C_Pi*Atom.LV[i].g*Q1*Nul*Mag.B*I;
                                }

                                if (InCol_Flag > 0 && K1==K2 &&Q1==Q2) {
                                    for (indx=0; indx<Atom.Num_Level; indx++) {
                                        if (i!=indx) {
                                            SEE_Matrix[a][b] -= InCol_Rates[i][indx] \
                                                *Density_Electron;
                                        }
                                    }
                                }
                            }else{
                                if (Atom.B[j][i]>0) {
                                    SEE_Matrix[a][b] = TE(i, Atom.LV[i].J, K1, Q1, j, \
    Atom.LV[j].J, K2, Q2, Atom)+TS(i, Atom.LV[i].J, K1, Q1, j, Atom.LV[j].J, K2, Q2, Atom);
                                }else{
                                    SEE_Matrix[a][b]=0;
                                }

                                if (InCol_Flag > 0 && K1==K2 &&Q1==Q2) {
                                    if (tmp2 != 0) {
                                        SEE_Matrix[a][b] += sqrt((2*Atom.LV[j].J+1) \
    /(2*Atom.LV[i].J+1))*tmp2*InCol_Rates[j][i]*Density_Electron;
                                    }
                                }
                            }
                            b++;
                        }
                    }
                }
                a++;
            }
        }
    }
    return;
}

static void SEE2(STRUCT_ATOM Atom, complex double **SEE_Matrix, double **InCol_Rates, \
                 int InCol_Flag, double Density_Electron, STRUCT_MAG Mag){
    
    /***********************************************************************************
     Purpose:
     Compute statistical equilibrium equations coefficients for multi-level atom.
     (only Rho with Q=0 and K even are taken into account)
     in the equation).
     Modified:
     1 Dec. 2019
     Input parameters:
     Atom, a structure saved the atom information.
     InCol_Rates, inelastic collisional rates.
     InCol_Flag, inelastic collision flag.
     Density_Electron, electron density.
     Mag, a structure saved the magnetic field.
     Output parameters:
     SEE_Matrix, the coefficient matrix.
     References:
     LL04, Chapter 7, Page 285, Eq 711.
     ***********************************************************************************/
    
    int i, j, K1, K2, Q1 = 0, Q2 = 0, indx, a = 0, b = 0;
    double tmp1 = 0, tmp2 = 0;
    
    for (i=0; i<Atom.Num_Level; i++) {
        for (K1=0; K1<=Atom.LV[i].J*2; K1+=2) {
            b = 0;
            for (j=0; j<Atom.Num_Level; j++) {
                if (InCol_Flag > 0) {
                    tmp1 = SJ(Atom.LV[i].J, Atom.LV[i].J, 0, Atom.LV[j].J, \
                              Atom.LV[j].J, 1);
                    
                    if (tmp1 != 0 ) {
                        tmp2 = pow(-1, K1)*SJ(Atom.LV[i].J, Atom.LV[i].J, K1, \
                                              Atom.LV[j].J, Atom.LV[j].J, 1)/
                        SJ(Atom.LV[i].J, Atom.LV[i].J, 0, Atom.LV[j].J, \
                           Atom.LV[j].J, 1);
                        
                    }else{
                        if (K1==0) {
                            tmp2 = 1;
                        }else{
                            tmp2 = 1;
                        }
                    }
                    
                }
                for (K2=0; K2<=Atom.LV[j].J*2; K2+=2) {
                    if (i>j) {
                        if (Atom.B[i][j]>0) {
                            SEE_Matrix[a][b] = TA(i, Atom.LV[i].J, K1, Q1, \
                                                  j, Atom.LV[j].J, K2, Q2, Atom);
                        }else{
                            SEE_Matrix[a][b]=0;
                        }
                        
                        if (InCol_Flag > 0 && K1==K2 && Q1==Q2) {
                            if (tmp2 != 0) {
                                SEE_Matrix[a][b] += sqrt((2*Atom.LV[j].J+1) \
    /(2*Atom.LV[i].J+1))*tmp2*InCol_Rates[j][i]*Density_Electron;
                            }
                        }
                    }else if(i==j){
                        SEE_Matrix[a][b] = -RA(i, Atom.LV[i].J, K1, Q1, K2, Q2, Atom) \
    -RE(i, Atom.LV[i].J, K1, Q1, K2, Q2, Atom)-RE(i, Atom.LV[i].J, K1, Q1, K2, Q2, Atom) \
    -RS(i, Atom.LV[i].J, K1, Q1, K2, Q2, Atom);
                        if (a==b && Q1!=0) {
                            SEE_Matrix[a][b] -= 2.*C_Pi*Atom.LV[i].g*Q1*Nul*Mag.B*I;
                        }
                        
                        if (InCol_Flag > 0 && K1==K2 &&Q1==Q2) {
                            for (indx=0; indx<Atom.Num_Level; indx++) {
                                if (i!=indx) {
                                    SEE_Matrix[a][b] -= InCol_Rates[i][indx] \
                                    *Density_Electron;
                                }
                            }
                        }
                    }else{
                        if (Atom.B[j][i]>0) {
                            SEE_Matrix[a][b] = TE(i, Atom.LV[i].J, K1, Q1, j, \
    Atom.LV[j].J, K2, Q2, Atom)+TS(i, Atom.LV[i].J, K1, Q1, j, Atom.LV[j].J, K2, Q2, Atom);
                        }else{
                            SEE_Matrix[a][b]=0;
                        }
                        
                        if (InCol_Flag > 0 && K1==K2 &&Q1==Q2) {
                            if (tmp2 != 0) {
                                SEE_Matrix[a][b] += sqrt((2*Atom.LV[j].J+1) \
    /(2*Atom.LV[i].J+1))*tmp2*InCol_Rates[j][i]*Density_Electron;
                            }
                        }
                    }
                    b++;
                }
            }
            a++;
        }
    }
    return;
}


void Rho_Compute(STRUCT_ATOM Atom, STRUCT_SPH_COORD Position, STRUCT_PARA Para, \
                 double **InCol_Rates, int InCol_Flag, int Flag_Reduced){
    
    /***********************************************************************************
     Purpose:
     Compute the multipolar components of the density matrix.
     Modified:
     1 Dec. 2019
     Input parameters:
     Atom, a structure saved the atom information.
     Position, position is the coronal model.
     Para, physical parameters at the position.
     InCol_Rates, inelastic collisional rates.
     InCol_Flag, inelastic collision flag.
     Flag_Reduced, if Flag_Cal, (only Rho with Q=0 and K even are taken into account
     in the equation).
     Output parameters:
     Atom, a structure saved the atom information.
     ***********************************************************************************/
    
    int i;
    complex double **Rho_Matrix, *EQ_B;
    Rho_Matrix = MATRIX_COMPLEX(0, Atom.Num_EQ-1, 0, Atom.Num_EQ-1, 1);
    EQ_B = VECTOR_COMPLEX(1, Atom.Num_EQ-1, 0);
    
    int *indx;
    indx = VECTOR_INT(1, Atom.Num_EQ-1, 1);
    
    if(Flag_Reduced){
        SEE2(Atom, Rho_Matrix, InCol_Rates, InCol_Flag, Para.Electron, Para.Mag);
    }else{
        SEE(Atom, Rho_Matrix, InCol_Rates, InCol_Flag, Para.Electron, Para.Mag);
    }

    for (i=1; i<Atom.Num_EQ; i++) {
        EQ_B[i] = -Rho_Matrix[i][0];
    }
  
    //do the LU decomposition
    ludcmp_complex(Rho_Matrix, Atom.Num_EQ-1, indx);
    
    //sove the equations output the solution in Rho_M
    lubksb_complex(Rho_Matrix, Atom.Num_EQ-1, indx, EQ_B);
    
    Atom.LV[0].Rho[0][0] = 1;
   
    int K, Q, ii = 1;
    
    if(Flag_Reduced){
        if (Atom.LV[0].J>=0) {
            for (K=1; K<=Atom.LV[0].J*2; K+=2) {
                Atom.LV[0].Rho[K][0] = EQ_B[ii];
                ii++;
            }
        }
        
        for (i=1; i<Atom.Num_Level; i++) {
            for (K=0; K<=Atom.LV[i].J*2; K+=2) {
                Atom.LV[i].Rho[K][0] = EQ_B[ii];
                ii++;
            }
        }
    }else{
        if (Atom.LV[0].J>=0) {
            for (K=1; K<=Atom.LV[0].J*2; K++) {
                for (Q=-K; Q<=K; Q++) {
                    Atom.LV[0].Rho[K][Q] = EQ_B[ii];
                    ii++;
                }
            }
        }
        
        for (i=1; i<Atom.Num_Level; i++) {
            for (K=0; K<=Atom.LV[i].J*2; K++) {
                for (Q=-K; Q<=K; Q++) {
                    Atom.LV[i].Rho[K][Q] = EQ_B[ii];
                    ii++;
                }
            }
        }
    }
   
    FREE_MATRIX_COMPLEX(Rho_Matrix, 0, 0);
    FREE_VECTOR_INT(indx, 1);
    FREE_VECTOR_COMPLEX(EQ_B, 1);
    return;
}


extern void Init_Output(STRUCT_OUTPUT Output, STRUCT_ATOM Atom){
    
    /***********************************************************************************
     Purpose:
     Initialize the output structure.
     Modified:
     1 Dec. 2019
     Input parameters:
     Atom, a structure saved the atom information.
     Output parameters:
     Output, the output structure.
     ***********************************************************************************/
    
    int i, j;
    double tmp;
    //LL04 Chaper 10 Eq.10.11
    for (i=0; i<Output.Num_transition; i++) {
        Output.TR[i].w[0] = 1;
        tmp = pow(-1, 1+Atom.LV[Output.TR[i].au].J+Atom.LV[Output.TR[i].al].J) \
            *sqrt(3.*(2.*Atom.LV[Output.TR[i].au].J+1));
        for (j=1; j<3; j++) {
            Output.TR[i].w[j] = tmp*SJ(1, 1, j, Atom.LV[Output.TR[i].au].J, \
                Atom.LV[Output.TR[i].au].J, Atom.LV[Output.TR[i].al].J);
        }
        Output.TR[i].Indx_Transition = -1;
        for (j=0; j<Atom.Num_Transition; j++) {
            if (Output.TR[i].au==Atom.TR[j].au && Output.TR[i].al==Atom.TR[j].al) {
                Output.TR[i].Indx_Transition = j;
                break;
            }
        }
        if (Output.TR[i].Indx_Transition < 0) {
            fprintf(stderr, "Can not find the input transiton \n");
            exit(1);
        }
    }
    return;
}


extern void Free_Output(STRUCT_OUTPUT *Output){
    
    /***********************************************************************************
     Purpose:
     Free the output structure.
     Modified:
     1 Dec. 2019
     Input parameters:
     Output, the output structure.
     ***********************************************************************************/
    
    free((*Output).TR);
    free(Output);
    return;
}


extern void EMCoefi_Inteall(STRUCT_ATOM Atom, STRUCT_PARA Para, STRUCT_OUTPUT output, \
                            int Flag_Cal, int Flag_Reduced){
    
    /***********************************************************************************
     Purpose:
     Compute the integrated emission coefficients for a transition (Q is along the
     z direction).
     Modified:
     1 Dec. 2019
     Input parameters:
     Atom, a structure saved the atom information.
     Para, physical parameters at the position.
     Output, the output structure.
     Flag_Cal, if Flag_Cal != 0 use Eq 10.31 in Chapter 10 Page 521, else use Eq 7.15 in
     Chapter 7, Page 289
     Flag_Reduced, if Flag_Cal, (only Rho with Q=0 and K even are taken into account
     in the equation).
     Output parameters:
     Output, the output structure.
     References:
     LL04
     ***********************************************************************************/
    
    int i, j, m;
    complex double ***T_KP, ***T_KQ;
    T_KP = MATRIX3_RHO(3, 2, 1);
    T_KQ = MATRIX3_RHO(3, 2, 1);
    
    TP_tensor(T_KP);
    
    int K, Q;
    for (i=0; i<4; i++) {
        for (K=1; K<=2; K++) {
            for (Q=-K; Q<=K; Q++) {
                T_KQ[i][K][Q] = 0;
            }
        }
        
        T_KQ[i][2][0] = 0;
        
        for (K=1; K<=2; K++) {
            for (Q=-K; Q<=K; Q++) {
                for (j = -2; j <= 2; j++) {
                    for (m = -2;  m<=2; m++) {
                        if(creal(T_KP[i][K][j])!=0 || cimag(T_KP[i][K][j])!=0){
                            T_KQ[i][K][Q] += T_KP[i][K][j] \
                                *Djmn(0, -0.5*C_Pi, 0.5*C_Pi, K, j, m) \
                                *Djmn(Para.Mag.PhiB, Para.Mag.ThetaB, 0, K, m, Q);
                        }
                    }
                }
            }
        }
        
    }
    
    T_KQ[0][0][0] = 1;
    
    for (i=0; i<output.Num_transition; i++) {
        EMCoefi_Inte(Atom, T_KQ, output.TR+i, Flag_Cal, Flag_Reduced);
    }
    
    FREE_MATRIX3_RHO(T_KP);
    FREE_MATRIX3_RHO(T_KQ);
    
    return;
}


extern void EMCoefi_Inte(STRUCT_ATOM Atom, complex double ***T_KQ, \
            STRUCT_TRANSITION_OUT *Transition, int Flag_Cal, int Flag_Reduced){
    
    /***********************************************************************************
     Purpose:
     Compute the integrated emission coefficients for a transition.
     Modified:
     1 Dec. 2019
     Input parameters:
     Atom, a structure saved the atom information.
     ***T_KQ, geometrical tensor.
     *Transition, a structure saved the transition information
     Flag_Cal, if Flag_Cal != 0 use Eq 10.31 in Chapter 10 Page 521, else use Eq 7.15 in
     Chapter 7, Page 289
     Flag_Reduced, if Flag_Cal, (only Rho with Q=0 and K even are taken into account
     in the equation).
     Output parameters:
     *Transition, emmission coefficients saved in the structure.
     References:
     LL04
     ***********************************************************************************/
    
    double nu = C_c*(Atom.LV[(*Transition).au].Energy-Atom.LV[(*Transition).al].Energy);
    double Ju = Atom.LV[(*Transition).au].J;
    int i, K, Q;
    double tmp;

    for (i=0; i<3; i++) {
        (*Transition).Emcoefi[i] = 0;
    }
    
    if (Flag_Cal==0) {
        tmp=C_h*nu/4/C_Pi*sqrt(2*Ju+1.)*Atom.A[(*Transition).au][(*Transition).al];
        if (Flag_Reduced) {
            for (i=0; i<3; i++) {
                for (K=0; K<=2; K+=2) {
                    (*Transition).Emcoefi[i] += tmp*(*Transition).w[K] \
                    *creal(T_KQ[i][K][0]*Atom.LV[(*Transition).au].Rho[K][0]);
                }
            }
        }else{
            for (i=0; i<3; i++) {
                for (K=0; K<=2; K+=2) {
                    for (Q=-K; Q<=K; Q++) {
                        (*Transition).Emcoefi[i] += tmp*(*Transition).w[K] \
                        *creal(T_KQ[i][K][Q]*Atom.LV[(*Transition).au].Rho[K][Q]);
                    }
                }
            }
        }

        for (i=0; i>3; i++) {
            (*Transition).Emcoefi[i] /= sqrt(2.*Atom.LV[0].J+1);
        }
        if ((*Transition).Flag_Mag_Diple>0) {
            for (i=1; i<3; i++) {
                (*Transition).Emcoefi[i] = -(*Transition).Emcoefi[i];
            }
        }
    }else{
        int Ku, Qu, q, qp, Mu2, Mup2, Ml2;
        double Mu, Mup, Ml, tmp1, tmp2;
        double Jl = Atom.LV[(*Transition).al].J;
        int Ju2 = (int)(2*Ju), Jl2=(int)(2*Jl);
        
        tmp = C_h*nu/4/C_Pi*(2*Ju+1.)*Atom.A[(*Transition).au][(*Transition).al];
    
        for (K=0; K<=2; K++) {
            for (Ku=0; Ku<=Ju2; Ku++) {
                tmp1 = tmp*sqrt(3.*(2*K+1.)*(2*Ku+1.));
                for (Mu2=-Ju2; Mu2<=Ju2; Mu2+=2) {
                    Mu = 0.5*Mu2;
                    for (Mup2=-Ju2; Mup2<=Ju2; Mup2+=2) {
                        Mup = 0.5* Mup2;
                        for (Ml2=-Jl2; Ml2<=Jl2; Ml2+=2) {
                            Ml = 0.5*Ml2;
                            for (q=-1; q<=1; q++) {
                                for (qp=-1; qp<=1; qp++) {
                                    for (Q=-K; Q<=K; Q++) {
                                        for (Qu=-Ku; Qu<=Ku; Qu++) {
                                            tmp2 = tmp1*pow(-1, 1+Ju-Mu+qp) \
        *TJ(Ju, Jl, 1, -Mu, Ml, -q)*TJ(Ju, Jl, 1, -Mup, Ml, -qp) \
        *TJ(1, 1, K, q, -qp, -Q)*TJ(Ju, Ju, Ku, Mup, -Mu, -Qu);
                                        
                                            for (i=0; i<3; i++) {
                                                (*Transition).Emcoefi[i] += \
        tmp2*creal(T_KQ[i][K][Q]*Atom.LV[(*Transition).au].Rho[Ku][Qu]);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        for (i=0; i>3; i++) {
            (*Transition).Emcoefi[i] /= sqrt(2.*Atom.LV[0].J+1);
        }
        if ((*Transition).Flag_Mag_Diple>0) {
            for (i=1; i<3; i++) {
                (*Transition).Emcoefi[i] = -(*Transition).Emcoefi[i];
            }
        }
    }
    return;
}


extern void Thom_Scat_van(double R, double *PB, double q){
    
    /***********************************************************************************
     Purpose:
     Compute the Thomson scattering (a and b coefficients) in van de Hulst 1950BAN.
     Modified:
     1 Dec. 2019
     Input parameters:
     R, Radius.
     q, limb darkening coefficient.
     Output parameters:
     PB[2], A and B in Eqs 11 and 12 (van de Hulst 1950BAN)
     References:
     van de Hulst 1950BAN
     ***********************************************************************************/
    
    double gamma = asin(1/R);
    double temp1, temp2;
    
    temp1 = (1-q)/(1-1/3.*q)*(2*(1-cos(gamma)))+q/(1-1/3.*q) \
        *(1-cos(gamma)*cos(gamma)/sin(gamma)*log((1+sin(gamma))/cos(gamma)));
    temp2 = (1-q)/(1-1/3.*q)*(2/3.*(1-cos(gamma)*cos(gamma)*cos(gamma)))+q/(1-1/3.*q) \
    *(1/4.+sin(gamma)*sin(gamma)/4.-cos(gamma)*cos(gamma)*cos(gamma)*cos(gamma)/4.\
      /sin(gamma)*log((1+sin(gamma))/cos(gamma)));
    
    PB[0] = (temp1+temp2)/4.;
    PB[1] = (temp1-temp2)/2.;
    
    return;
}


extern void Thom_Scat(double r, double *PB, double u1, double u2){
    
    /***********************************************************************************
     Purpose:
     Compute the Thomson scattering (a and b coefficients) in van de Hulst 1950BAN.
     Modified:
     1 Dec. 2019
     Input parameters:
     R, Radius.
     u1, u2, limb darkening coefficient.
     Output parameters:
     PB[2], a and b in Eqs 11 and 12 (van de Hulst 1950BAN).
     References:
     LL04
     Caution:
     a and b are difined in van de Hulst 1950 a factor of pi is not included.
     ***********************************************************************************/

    double  Sr, Cr, a[3], b[3], Jv, Kv;
    
    Sr = 1/r;
    Cr = sqrt(1-Sr*Sr);
    a[0] = 1-Cr;
    a[1] = Cr-0.5-0.5*Cr*Cr/Sr*log((1+Sr)/Cr);
    a[2] = (Cr+2)*(Cr-1)/3/(Cr+1);
    
    b[0] = 1./3.*(1-Cr*Cr*Cr);
    b[1] = 1./24.*(8*Cr*Cr*Cr-3*Cr*Cr-2)-1./8.*Cr*Cr*Cr*Cr/Sr*log((1+Sr)/Cr);
    b[2] = (Cr-1)*(3*Cr*Cr*Cr+6*Cr*Cr+4*Cr+2)/15./(Cr+1);
    
    Jv = 0.5*(a[0]+a[1]*u1+a[2]*u2);
    Kv = 0.5*(b[0]+b[1]*u1+b[2]*u2);
        
    PB[0] = (Jv+Kv);
    PB[1] = (Jv-Kv)*2.;
    
    return;
}

/**********************************************************************************************/
