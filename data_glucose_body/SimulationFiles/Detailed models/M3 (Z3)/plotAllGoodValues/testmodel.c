/* General includes */
#include "math.h"
#include "stddef.h"
#include "stdarg.h"

/* Integrator interface include */
#include "interfaceSER.h"

/* MODEL RELATED include */
#include "testmodel.h"

/* SBTOOLBOX relevant includes */
#include "mexODEDATA.h"

void model(double time, double *statevector, double *ODERHSvector, ModelData *modeldataPtr, int CVODEflag, double *variablevector, double *reactionvector, double *gout, int *eventvector)
{
    double r0,r1,r2,r11,r12,r22,r1x2,r11x2,r1x22,r1x22d,r11x22,rend,rendP,iendIR,iend,rPbasal,IRS,IRSiP,X,X_P,PI3K,PI3K_,PDK1,PDK1_,PKC,PKC_P,PKB,PKB_P,mTOR,mTOR_,GLUT4_C,GLUT4_M;
    double ins,glucose,a1,a2,d1,d2,Kcr,Kex,Kend,Kdp,Kcat,Km,kfbasal,krbasal,k21,k22,k23,k24,k2b,k3f,k3b,k4f,k4b,k5f,k5b,k6f,k6b,k7f,k7b,k8f,k8b,k91,k92,k9b,k5Basic,k5BasicWb,k_glut4,k_glut1,KmG1,KmG4,kbf;
    double MeasIR,MeasIRS,MeasPKB,MeasAllIntIR,totX,itot,KD,S2,S1,K4,K8;
    double vglucoseuptake,vglut1,vglut4,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15,R16,R17,R18,R19,R20,R21,R22,R23,R24,R25,R26,R27,R28,R29,R30,R31,R32,R33,R34,R35,R36,R37,R38,R39,R40,R41,R42,R43,R44,R45,R46,R47;
    double R48,v2f,v2b,v3f,v3b,v4f,v4b,v5f,v5b,v6f,v6b,v7f,v7b,v8f,v8b,v9f,v9b;

    /* STATES */
    r0 = statevector[0];
    r1 = statevector[1];
    r2 = statevector[2];
    r11 = statevector[3];
    r12 = statevector[4];
    r22 = statevector[5];
    r1x2 = statevector[6];
    r11x2 = statevector[7];
    r1x22 = statevector[8];
    r1x22d = statevector[9];
    r11x22 = statevector[10];
    rend = statevector[11];
    rendP = statevector[12];
    iendIR = statevector[13];
    iend = statevector[14];
    rPbasal = statevector[15];
    IRS = statevector[16];
    IRSiP = statevector[17];
    X = statevector[18];
    X_P = statevector[19];
    PI3K = statevector[20];
    PI3K_ = statevector[21];
    PDK1 = statevector[22];
    PDK1_ = statevector[23];
    PKC = statevector[24];
    PKC_P = statevector[25];
    PKB = statevector[26];
    PKB_P = statevector[27];
    mTOR = statevector[28];
    mTOR_ = statevector[29];
    GLUT4_C = statevector[30];
    GLUT4_M = statevector[31];
    /* PARAMETERS */
    ins = modeldataPtr->parametervector[0]; /* 0 */
    glucose = modeldataPtr->parametervector[1]; /* 0 */
    a1 = modeldataPtr->parametervector[2]; /* 4.48251e+008 */
    a2 = modeldataPtr->parametervector[3]; /* 4.32189e+006 */
    d1 = modeldataPtr->parametervector[4]; /* 0.772261 */
    d2 = modeldataPtr->parametervector[5]; /* 0.0122058 */
    Kcr = modeldataPtr->parametervector[6]; /* 0.00136404 */
    Kex = modeldataPtr->parametervector[7]; /* 37.0819 */
    Kend = modeldataPtr->parametervector[8]; /* 30.6825 */
    Kdp = modeldataPtr->parametervector[9]; /* 0.000950083 */
    Kcat = modeldataPtr->parametervector[10]; /* 237.519 */
    Km = modeldataPtr->parametervector[11]; /* 3.01819 */
    kfbasal = modeldataPtr->parametervector[12]; /* 0.497522 */
    krbasal = modeldataPtr->parametervector[13]; /* 128043 */
    k21 = modeldataPtr->parametervector[14]; /* 0.00964586 */
    k22 = modeldataPtr->parametervector[15]; /* 2374.98 */
    k23 = modeldataPtr->parametervector[16]; /* 8.92e-008 */
    k24 = modeldataPtr->parametervector[17]; /* 1.168e-007 */
    k2b = modeldataPtr->parametervector[18]; /* 608.584 */
    k3f = modeldataPtr->parametervector[19]; /* 8.11194 */
    k3b = modeldataPtr->parametervector[20]; /* 0.18953 */
    k4f = modeldataPtr->parametervector[21]; /* 180474 */
    k4b = modeldataPtr->parametervector[22]; /* 1.0309e+006 */
    k5f = modeldataPtr->parametervector[23]; /* 376367 */
    k5b = modeldataPtr->parametervector[24]; /* 5625.18 */
    k6f = modeldataPtr->parametervector[25]; /* 1.04646e+007 */
    k6b = modeldataPtr->parametervector[26]; /* 2380 */
    k7f = modeldataPtr->parametervector[27]; /* 1859.07 */
    k7b = modeldataPtr->parametervector[28]; /* 150488 */
    k8f = modeldataPtr->parametervector[29]; /* 3.34849e+006 */
    k8b = modeldataPtr->parametervector[30]; /* 35149.7 */
    k91 = modeldataPtr->parametervector[31]; /* 1.45e-008 */
    k92 = modeldataPtr->parametervector[32]; /* 0.00710824 */
    k9b = modeldataPtr->parametervector[33]; /* 3.96307 */
    k5Basic = modeldataPtr->parametervector[34]; /* 0.00349967 */
    k5BasicWb = modeldataPtr->parametervector[35]; /* 1.19e-008 */
    k_glut4 = modeldataPtr->parametervector[36]; /* 59.1318 */
    k_glut1 = modeldataPtr->parametervector[37]; /* 0.0816969 */
    KmG1 = modeldataPtr->parametervector[38]; /* 157.063 */
    KmG4 = modeldataPtr->parametervector[39]; /* 71.6209 */
    kbf = modeldataPtr->parametervector[40]; /* 0.00375236 */
    /* VARIABLES */
    MeasIR = r1x2+r11x2+r1x22+r1x22d+r11x22+rendP+rPbasal;
    MeasIRS = IRSiP;
    MeasPKB = PKB_P;
    MeasAllIntIR = rend+rendP;
    totX = X+X_P;
    itot = ins*1e-9;
    KD = 7e-6;
    S2 = (4.0*itot+KD-sqrt(pow(KD,2.0)+8.0*itot*KD))/8.0;
    S1 = itot-S2;
    K4 = 1400.0;
    K8 = 0.01;
    /* REACTIONS */
    vglucoseuptake = k_glut1*glucose/(KmG1+glucose)+k_glut4*GLUT4_M*glucose/(KmG4+glucose);
    vglut1 = k_glut1*glucose/(KmG1+glucose);
    vglut4 = k_glut4*GLUT4_M*glucose/(KmG4+glucose);
    R1 = 2.0*a1*S1*r0;
    R2 = 2.0*a2*S1*r0;
    R3 = a1*S1*r1;
    R4 = a1*S1*r2;
    R5 = d1*r1;
    R6 = a2*S1*r1;
    R7 = a2*S1*r2;
    R8 = d2*r2;
    R9 = Kcr*r1;
    R10 = Kcr*r2;
    R11 = a1*S1*r1x2;
    R12 = 2.0*d1*r11;
    R13 = d1*r12;
    R14 = a2*S1*r1x2;
    R15 = d2*r12;
    R16 = 2.0*d2*r22;
    R17 = 2.0*Kcr*r11;
    R18 = Kcr*r12;
    R19 = d2*r1x2;
    R20 = Kcr*r12;
    R21 = 2.0*Kcr*r22;
    R22 = d1*r1x2;
    R23 = a2*S2*r1x2;
    R24 = d1*r11x2;
    R25 = d2*r1x22;
    R26 = d2*r11x2;
    R27 = d2*r1x22;
    R28 = d1*r11x2;
    R29 = d1*r1x22;
    R30 = a1*S1*r1x22;
    R31 = a2*S1*r11x2;
    R32 = K4*S1*r1x22;
    R33 = K8*r1x22d;
    R34 = d2*r1x22d;
    R35 = d1*r11x22;
    R36 = d2*r11x22;
    R37 = Kex*rend;
    R38 = Kex*iend;
    R39 = (Kend)*r1x2;
    R40 = (Kend)*r11x2;
    R41 = (Kend)*r1x22;
    R42 = (Kend)*r1x22d;
    R43 = (Kend)*r11x22;
    R44 = (Kdp+Kcat*(X_P/totX)/(Km+(X_P/totX)))*rendP;
    R45 = (Kdp+Kcat*(X_P/totX)/(Km+(X_P/totX)))*iendIR;
    R46 = kfbasal*r0;
    R47 = krbasal*rPbasal;
    R48 = Kend*rPbasal;
    v2f = k21*IRS*((r1x2+r11x2+r1x22+r1x22d+r11x22+rPbasal)+k22*rendP)*(1.0+k23*PKC_P+k24*mTOR);
    v2b = k2b*IRSiP;
    v3f = k3f*X*IRSiP;
    v3b = k3b*X_P;
    v4f = k4f*PI3K*IRSiP;
    v4b = k4b*PI3K_;
    v5f = k5f*PDK1*PI3K_;
    v5b = k5b*PDK1_;
    v6f = k6f*PKC*PDK1_;
    v6b = k6b*PKC_P;
    v7f = k7f*PKB*PDK1_;
    v7b = k7b*PKB_P;
    v8f = k8f*mTOR*PKB_P;
    v8b = k8b*mTOR_;
    v9f = k91*GLUT4_C*PKC_P+k92*GLUT4_C*PKB_P+k5BasicWb*GLUT4_C;
    v9b = k9b*GLUT4_M;
    /* Dependent on the CVODEflag, either the ODERHSvector are determined, */
    /* or the reaction and variable values. */
    if (CVODEflag == CVODE_RHS) {
        /* DIFFERENTIAL EQUATIONS (OUTPUT) */
        ODERHSvector[0] = -R1-R2+R5+R8+R37-R46+R47;
        ODERHSvector[1] = +R1-R3-R5-R6-R9+R12+R15+R19;
        ODERHSvector[2] = +R2-R4-R7-R8-R10+R13+R16+R22;
        ODERHSvector[3] = +R3-R12-R17+R26;
        ODERHSvector[4] = +R4+R6-R13-R15-R18-R20+R27+R28;
        ODERHSvector[5] = +R7-R16-R21+R29;
        ODERHSvector[6] = +R9+R10-R11-R14-R19-R22-R23+R24+R25+R34-R39;
        ODERHSvector[7] = +R11+R17+R20-R24-R26-R28-R31+R36-R40;
        ODERHSvector[8] = +R14+R18+R21-R25-R27-R29-R30-R32+R33+R35-R41;
        ODERHSvector[9] = +R23+R32-R33-R34-R42;
        ODERHSvector[10] = +R30+R31-R35-R36-R43;
        ODERHSvector[11] = -R37+R44;
        ODERHSvector[12] = -R44+R39+R40+R41+R42+R43+R48;
        ODERHSvector[13] = +R39+2.0*R40+2.0*R41+3.0*R42+3.0*R43-R45;
        ODERHSvector[14] = -R38+R45;
        ODERHSvector[15] = R46-R47-R48;
        ODERHSvector[16] = v2b-v2f;
        ODERHSvector[17] = -v2b+v2f;
        ODERHSvector[18] = v3b-v3f;
        ODERHSvector[19] = -v3b+v3f;
        ODERHSvector[20] = v4b-v4f;
        ODERHSvector[21] = -v4b+v4f;
        ODERHSvector[22] = v5b-v5f;
        ODERHSvector[23] = -v5b+v5f;
        ODERHSvector[24] = v6b-v6f;
        ODERHSvector[25] = -v6b+v6f;
        ODERHSvector[26] = v7b-v7f;
        ODERHSvector[27] = -v7b+v7f;
        ODERHSvector[28] = v8b-v8f;
        ODERHSvector[29] = -v8b+v8f;
        ODERHSvector[30] = v9b-v9f;
        ODERHSvector[31] = -v9b+v9f;
    } else if (CVODEflag == CVODE_VARREAC) {
        variablevector[0] = MeasIR;
        variablevector[1] = MeasIRS;
        variablevector[2] = MeasPKB;
        variablevector[3] = MeasAllIntIR;
        variablevector[4] = totX;
        variablevector[5] = itot;
        variablevector[6] = KD;
        variablevector[7] = S2;
        variablevector[8] = S1;
        variablevector[9] = K4;
        variablevector[10] = K8;
        reactionvector[0] = vglucoseuptake;
        reactionvector[1] = vglut1;
        reactionvector[2] = vglut4;
        reactionvector[3] = R1;
        reactionvector[4] = R2;
        reactionvector[5] = R3;
        reactionvector[6] = R4;
        reactionvector[7] = R5;
        reactionvector[8] = R6;
        reactionvector[9] = R7;
        reactionvector[10] = R8;
        reactionvector[11] = R9;
        reactionvector[12] = R10;
        reactionvector[13] = R11;
        reactionvector[14] = R12;
        reactionvector[15] = R13;
        reactionvector[16] = R14;
        reactionvector[17] = R15;
        reactionvector[18] = R16;
        reactionvector[19] = R17;
        reactionvector[20] = R18;
        reactionvector[21] = R19;
        reactionvector[22] = R20;
        reactionvector[23] = R21;
        reactionvector[24] = R22;
        reactionvector[25] = R23;
        reactionvector[26] = R24;
        reactionvector[27] = R25;
        reactionvector[28] = R26;
        reactionvector[29] = R27;
        reactionvector[30] = R28;
        reactionvector[31] = R29;
        reactionvector[32] = R30;
        reactionvector[33] = R31;
        reactionvector[34] = R32;
        reactionvector[35] = R33;
        reactionvector[36] = R34;
        reactionvector[37] = R35;
        reactionvector[38] = R36;
        reactionvector[39] = R37;
        reactionvector[40] = R38;
        reactionvector[41] = R39;
        reactionvector[42] = R40;
        reactionvector[43] = R41;
        reactionvector[44] = R42;
        reactionvector[45] = R43;
        reactionvector[46] = R44;
        reactionvector[47] = R45;
        reactionvector[48] = R46;
        reactionvector[49] = R47;
        reactionvector[50] = R48;
        reactionvector[51] = v2f;
        reactionvector[52] = v2b;
        reactionvector[53] = v3f;
        reactionvector[54] = v3b;
        reactionvector[55] = v4f;
        reactionvector[56] = v4b;
        reactionvector[57] = v5f;
        reactionvector[58] = v5b;
        reactionvector[59] = v6f;
        reactionvector[60] = v6b;
        reactionvector[61] = v7f;
        reactionvector[62] = v7b;
        reactionvector[63] = v8f;
        reactionvector[64] = v8b;
        reactionvector[65] = v9f;
        reactionvector[66] = v9b;
    } else if (CVODEflag == CVODE_EVENTS) {
        /* EVENTS */
    } else if (CVODEflag == CVODE_EVENTASSIGNMENT) {
        /* EVENT ASSIGNMENTS */
    }
}
