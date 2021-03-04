/* General includes */
#include "math.h"
#include "stddef.h"
#include "stdarg.h"

/* Integrator interface include */
#include "interfaceSER.h"

/* MODEL RELATED include */
#include "testmodelWb.h"

/* SBTOOLBOX relevant includes */
#include "mexODEDATA.h"

void model(double time, double *statevector, double *ODERHSvector, ModelData *modeldataPtr, int CVODEflag, double *variablevector, double *reactionvector, double *gout, int *eventvector)
{
    double IR,IRins,IR_P,IRiP,IRi,IRS,IRSiP,X,X_P,PKB,PKB_P,GLUT4_C,GLUT4_M,G_p,G_t,I_l,I_p,Q_sto1,Q_sto2,Q_gut,I_1,I_d,INS,I_po,Y;
    double ins,glucose,k1a,k1aBasic,k1b,k1c,k1d,k1e,k1f,k1g,k1r,k21,k22,k2b,k3f,k3b,k4f,k4b,k5f,k5b,k_glut4,k_glut1,KmG1,KmG4,V_G,k_1,k_2,G_b,V_I,m_1,m_2,m_4,m_5,m_6,HE_b,I_b,S_b,S_b_minus,k_max,k_min,k_abs,k_gri,f,b,d,BW,k_p1,k_p2,k_p3,k_p4;
    double k_i,U_ii,V_m0,V_mX,K_m0,V_f0,V_fX,K_f0,p_2U,part,K,alpha,beta,gamma,k_e1,k_e2,D;
    double MeasIR,MeasIRS,MeasPKB,aa,cc,EGP,V_mmax,V_fmax,E,S,U_idf,I,G,HE,m_3,Q_sto,Ra,k_empt,U_idm,U_id,U,S_po;
    double glucoseuptake,glut1,glut4,v1a,v1b,v1c,v1d,v1e,v1g,v1r,v2f,v2b,v3f,v3b,v4f,v4b,v5f,v5b;

    /* STATES */
    IR = statevector[0];
    IRins = statevector[1];
    IR_P = statevector[2];
    IRiP = statevector[3];
    IRi = statevector[4];
    IRS = statevector[5];
    IRSiP = statevector[6];
    X = statevector[7];
    X_P = statevector[8];
    PKB = statevector[9];
    PKB_P = statevector[10];
    GLUT4_C = statevector[11];
    GLUT4_M = statevector[12];
    G_p = statevector[13];
    G_t = statevector[14];
    I_l = statevector[15];
    I_p = statevector[16];
    Q_sto1 = statevector[17];
    Q_sto2 = statevector[18];
    Q_gut = statevector[19];
    I_1 = statevector[20];
    I_d = statevector[21];
    INS = statevector[22];
    I_po = statevector[23];
    Y = statevector[24];
    /* PARAMETERS */
    ins = modeldataPtr->parametervector[0]; /* 0 */
    glucose = modeldataPtr->parametervector[1]; /* 17 */
    k1a = modeldataPtr->parametervector[2]; /* 0.0675116 */
    k1aBasic = modeldataPtr->parametervector[3]; /* 0.00920433 */
    k1b = modeldataPtr->parametervector[4]; /* 69.2098 */
    k1c = modeldataPtr->parametervector[5]; /* 0.340865 */
    k1d = modeldataPtr->parametervector[6]; /* 0.00368488 */
    k1e = modeldataPtr->parametervector[7]; /* 0.118251 */
    k1f = modeldataPtr->parametervector[8]; /* 0.150493 */
    k1g = modeldataPtr->parametervector[9]; /* 1896.15 */
    k1r = modeldataPtr->parametervector[10]; /* 4857.43 */
    k21 = modeldataPtr->parametervector[11]; /* 858.671 */
    k22 = modeldataPtr->parametervector[12]; /* 31102.7 */
    k2b = modeldataPtr->parametervector[13]; /* 6351.66 */
    k3f = modeldataPtr->parametervector[14]; /* 153767 */
    k3b = modeldataPtr->parametervector[15]; /* 18.197 */
    k4f = modeldataPtr->parametervector[16]; /* 1798.87 */
    k4b = modeldataPtr->parametervector[17]; /* 928.917 */
    k5f = modeldataPtr->parametervector[18]; /* 0.331575 */
    k5b = modeldataPtr->parametervector[19]; /* 2.1637 */
    k_glut4 = modeldataPtr->parametervector[20]; /* 8.91135 */
    k_glut1 = modeldataPtr->parametervector[21]; /* 0.112454 */
    KmG1 = modeldataPtr->parametervector[22]; /* 1.00417 */
    KmG4 = modeldataPtr->parametervector[23]; /* 5100.22 */
    V_G = modeldataPtr->parametervector[24]; /* 1.88 */
    k_1 = modeldataPtr->parametervector[25]; /* 0.065 */
    k_2 = modeldataPtr->parametervector[26]; /* 0.079 */
    G_b = modeldataPtr->parametervector[27]; /* 95 */
    V_I = modeldataPtr->parametervector[28]; /* 0.05 */
    m_1 = modeldataPtr->parametervector[29]; /* 0.19 */
    m_2 = modeldataPtr->parametervector[30]; /* 0.484 */
    m_4 = modeldataPtr->parametervector[31]; /* 0.194 */
    m_5 = modeldataPtr->parametervector[32]; /* 0.0304 */
    m_6 = modeldataPtr->parametervector[33]; /* 0.6471 */
    HE_b = modeldataPtr->parametervector[34]; /* 0.6 */
    I_b = modeldataPtr->parametervector[35]; /* 25 */
    S_b = modeldataPtr->parametervector[36]; /* 1.8 */
    S_b_minus = modeldataPtr->parametervector[37]; /* -1.8 */
    k_max = modeldataPtr->parametervector[38]; /* 0.0558 */
    k_min = modeldataPtr->parametervector[39]; /* 0.008 */
    k_abs = modeldataPtr->parametervector[40]; /* 0.057 */
    k_gri = modeldataPtr->parametervector[41]; /* 0.0558 */
    f = modeldataPtr->parametervector[42]; /* 0.9 */
    b = modeldataPtr->parametervector[43]; /* 0.82 */
    d = modeldataPtr->parametervector[44]; /* 0.01 */
    BW = modeldataPtr->parametervector[45]; /* 78 */
    k_p1 = modeldataPtr->parametervector[46]; /* 2.7 */
    k_p2 = modeldataPtr->parametervector[47]; /* 0.0021 */
    k_p3 = modeldataPtr->parametervector[48]; /* 0.009 */
    k_p4 = modeldataPtr->parametervector[49]; /* 0.0618 */
    k_i = modeldataPtr->parametervector[50]; /* 0.0079 */
    U_ii = modeldataPtr->parametervector[51]; /* 1 */
    V_m0 = modeldataPtr->parametervector[52]; /* 2.5 */
    V_mX = modeldataPtr->parametervector[53]; /* 0.047 */
    K_m0 = modeldataPtr->parametervector[54]; /* 225.59 */
    V_f0 = modeldataPtr->parametervector[55]; /* 2.5 */
    V_fX = modeldataPtr->parametervector[56]; /* 0.047 */
    K_f0 = modeldataPtr->parametervector[57]; /* 225.59 */
    p_2U = modeldataPtr->parametervector[58]; /* 0.0331 */
    part = modeldataPtr->parametervector[59]; /* 0.2 */
    K = modeldataPtr->parametervector[60]; /* 2.3 */
    alpha = modeldataPtr->parametervector[61]; /* 0.05 */
    beta = modeldataPtr->parametervector[62]; /* 0.11 */
    gamma = modeldataPtr->parametervector[63]; /* 0.5 */
    k_e1 = modeldataPtr->parametervector[64]; /* 0.0005 */
    k_e2 = modeldataPtr->parametervector[65]; /* 339 */
    D = modeldataPtr->parametervector[66]; /* 78000 */
    /* VARIABLES */
    MeasIR = IR_P+IRiP;
    MeasIRS = IRSiP;
    MeasPKB = PKB_P;
    aa = 5.0/2.0/(1.0-b)/D;
    cc = 5.0/2.0/d/D;
    EGP = k_p1-k_p2*G_p-k_p3*I_d-k_p4*I_po;
    V_mmax = (1.0-part)*(V_m0+V_mX*INS);
    V_fmax = part*(V_f0+V_fX*INS);
    E = 0.0;
    S = gamma*I_po;
    U_idf = V_fmax*G_t/(K_f0+G_t);
    I = I_p/V_I;
    G = G_p/V_G;
    HE = (-m_5*S)+m_6;
    m_3 = HE*m_1/(1.0-HE);
    Q_sto = Q_sto1+Q_sto2;
    Ra = f*k_abs*Q_gut/BW;
    k_empt = k_min+(k_max-k_min)/2.0*(tanh(aa*(Q_sto-b*D))-tanh(cc*(Q_sto-d*D))+2.0);
    U_idm = V_mmax*G_t/(K_m0+G_t);
    U_id = U_idm+U_idf;
    U = U_ii+U_id;
    S_po = Y+K*(EGP+Ra-E-U_ii-k_1*G_p+k_2*G_t)/V_G+S_b;
    /* REACTIONS */
    glucoseuptake = k_glut1*G_t/(KmG1+G_t)+k_glut4*GLUT4_M*G_t/(KmG4+G_t);
    glut1 = k_glut1*G_t/(KmG1+G_t);
    glut4 = k_glut4*GLUT4_M*G_t/(KmG4+G_t);
    v1a = k1a*(INS+5.0)*IR+k1aBasic*IR;
    v1b = k1b*IRins;
    v1c = k1c*IRins;
    v1d = k1d*IR_P;
    v1e = IRiP*(k1e+k1f*X_P/(1.0+X_P));
    v1g = k1g*IR_P;
    v1r = k1r*IRi;
    v2f = k21*(IR_P+k22*IRiP)*IRS;
    v2b = k2b*IRSiP;
    v3f = k3f*IRSiP*X;
    v3b = k3b*X_P;
    v4f = k4f*IRSiP*PKB;
    v4b = k4b*PKB_P;
    v5f = k5f*PKB_P*GLUT4_C;
    v5b = k5b*GLUT4_M;
    /* Dependent on the CVODEflag, either the ODERHSvector are determined, */
    /* or the reaction and variable values. */
    if (CVODEflag == CVODE_RHS) {
        /* DIFFERENTIAL EQUATIONS (OUTPUT) */
        ODERHSvector[0] = -v1a+v1b+v1r+v1g;
        ODERHSvector[1] = v1a-v1b-v1c;
        ODERHSvector[2] = v1c-v1d-v1g;
        ODERHSvector[3] = v1d-v1e;
        ODERHSvector[4] = v1e-v1r;
        ODERHSvector[5] = v2b-v2f;
        ODERHSvector[6] = -v2b+v2f;
        ODERHSvector[7] = v3b-v3f;
        ODERHSvector[8] = -v3b+v3f;
        ODERHSvector[9] = v4b-v4f;
        ODERHSvector[10] = -v4b+v4f;
        ODERHSvector[11] = v5b-v5f;
        ODERHSvector[12] = -v5b+v5f;
        ODERHSvector[13] = EGP+Ra-E-U_ii-k_1*G_p+k_2*G_t;
        ODERHSvector[14] = (-U_id)+k_1*G_p-k_2*G_t;
        ODERHSvector[15] = (-m_1*I_l)-m_3*I_l+m_2*I_p+S;
        ODERHSvector[16] = (-m_2*I_p)-m_4*I_p+m_1*I_l;
        ODERHSvector[17] = -k_gri*Q_sto1;
        ODERHSvector[18] = (-k_empt*Q_sto2)+k_gri*Q_sto1;
        ODERHSvector[19] = (-k_abs*Q_gut)+k_empt*Q_sto2;
        ODERHSvector[20] = -k_i*(I_1-I);
        ODERHSvector[21] = -k_i*(I_d-I_1);
        ODERHSvector[22] = (-p_2U*INS)+p_2U*(I-I_b);
        ODERHSvector[23] = (-gamma*I_po)+S_po;
        ODERHSvector[24] = -alpha*(Y-beta*(G-G_b));
    } else if (CVODEflag == CVODE_VARREAC) {
        variablevector[0] = MeasIR;
        variablevector[1] = MeasIRS;
        variablevector[2] = MeasPKB;
        variablevector[3] = aa;
        variablevector[4] = cc;
        variablevector[5] = EGP;
        variablevector[6] = V_mmax;
        variablevector[7] = V_fmax;
        variablevector[8] = E;
        variablevector[9] = S;
        variablevector[10] = U_idf;
        variablevector[11] = I;
        variablevector[12] = G;
        variablevector[13] = HE;
        variablevector[14] = m_3;
        variablevector[15] = Q_sto;
        variablevector[16] = Ra;
        variablevector[17] = k_empt;
        variablevector[18] = U_idm;
        variablevector[19] = U_id;
        variablevector[20] = U;
        variablevector[21] = S_po;
        reactionvector[0] = glucoseuptake;
        reactionvector[1] = glut1;
        reactionvector[2] = glut4;
        reactionvector[3] = v1a;
        reactionvector[4] = v1b;
        reactionvector[5] = v1c;
        reactionvector[6] = v1d;
        reactionvector[7] = v1e;
        reactionvector[8] = v1g;
        reactionvector[9] = v1r;
        reactionvector[10] = v2f;
        reactionvector[11] = v2b;
        reactionvector[12] = v3f;
        reactionvector[13] = v3b;
        reactionvector[14] = v4f;
        reactionvector[15] = v4b;
        reactionvector[16] = v5f;
        reactionvector[17] = v5b;
    } else if (CVODEflag == CVODE_EVENTS) {
        /* EVENTS */
    } else if (CVODEflag == CVODE_EVENTASSIGNMENT) {
        /* EVENT ASSIGNMENTS */
    }
}
