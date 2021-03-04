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
    double IR,IR_P,IRS,IRS_P,PKB,PKB_P,GLUT4_C,GLUT4_M,G_p,G_t,I_l,I_p,Q_sto1,Q_sto2,Q_gut,I_1,I_d,INS,I_po,Y;
    double ins,glucose,k1aBasic,k1f,k1b,k2f,k2b,k4f,k4b,k5Basic,k5BasicWb,k5f,k5b,k_glut4,k_glut4Wb,k_glut1,KmG1,KmG4,V_G,k_1,k_2,G_b,V_I,m_1,m_2,m_4,m_5,m_6,HE_b,I_b,S_b,S_b_minus,k_max,k_min,k_abs,k_gri,f,b,d,BW,k_p1,k_p2,k_p3,k_p4,k_i,U_ii,V_m0,V_mX,K_m0,V_f0;
    double V_fX,K_f0,p_2U,part,K,alpha,beta,gamma,k_e1,k_e2,D;
    double MeasIR,MeasIRS,MeasPKB,aa,cc,EGP,V_mmax,V_fmax,E,S,U_idf,I,G,HE,m_3,Q_sto,Ra,k_empt,U_idm,U_id,U,S_po;
    double glucoseuptake,glut1,glut4,v1f,v1b,v2f,v2b,v4f,v4b,v5f,v5b;

    /* STATES */
    IR = statevector[0];
    IR_P = statevector[1];
    IRS = statevector[2];
    IRS_P = statevector[3];
    PKB = statevector[4];
    PKB_P = statevector[5];
    GLUT4_C = statevector[6];
    GLUT4_M = statevector[7];
    G_p = statevector[8];
    G_t = statevector[9];
    I_l = statevector[10];
    I_p = statevector[11];
    Q_sto1 = statevector[12];
    Q_sto2 = statevector[13];
    Q_gut = statevector[14];
    I_1 = statevector[15];
    I_d = statevector[16];
    INS = statevector[17];
    I_po = statevector[18];
    Y = statevector[19];
    /* PARAMETERS */
    ins = modeldataPtr->parametervector[0]; /* 0 */
    glucose = modeldataPtr->parametervector[1]; /* 17 */
    k1aBasic = modeldataPtr->parametervector[2]; /* 107610 */
    k1f = modeldataPtr->parametervector[3]; /* 720.746 */
    k1b = modeldataPtr->parametervector[4]; /* 726898 */
    k2f = modeldataPtr->parametervector[5]; /* 479.078 */
    k2b = modeldataPtr->parametervector[6]; /* 999792 */
    k4f = modeldataPtr->parametervector[7]; /* 119471 */
    k4b = modeldataPtr->parametervector[8]; /* 159236 */
    k5Basic = modeldataPtr->parametervector[9]; /* 0.001 */
    k5BasicWb = modeldataPtr->parametervector[10]; /* 0.001 */
    k5f = modeldataPtr->parametervector[11]; /* 30433.9 */
    k5b = modeldataPtr->parametervector[12]; /* 47736.8 */
    k_glut4 = modeldataPtr->parametervector[13]; /* 79749.3 */
    k_glut4Wb = modeldataPtr->parametervector[14]; /* 79749.3 */
    k_glut1 = modeldataPtr->parametervector[15]; /* 1.6 */
    KmG1 = modeldataPtr->parametervector[16]; /* 294 */
    KmG4 = modeldataPtr->parametervector[17]; /* 21 */
    V_G = modeldataPtr->parametervector[18]; /* 1.88 */
    k_1 = modeldataPtr->parametervector[19]; /* 0.065 */
    k_2 = modeldataPtr->parametervector[20]; /* 0.079 */
    G_b = modeldataPtr->parametervector[21]; /* 95 */
    V_I = modeldataPtr->parametervector[22]; /* 0.05 */
    m_1 = modeldataPtr->parametervector[23]; /* 0.19 */
    m_2 = modeldataPtr->parametervector[24]; /* 0.484 */
    m_4 = modeldataPtr->parametervector[25]; /* 0.194 */
    m_5 = modeldataPtr->parametervector[26]; /* 0.0304 */
    m_6 = modeldataPtr->parametervector[27]; /* 0.6471 */
    HE_b = modeldataPtr->parametervector[28]; /* 0.6 */
    I_b = modeldataPtr->parametervector[29]; /* 25 */
    S_b = modeldataPtr->parametervector[30]; /* 1.8 */
    S_b_minus = modeldataPtr->parametervector[31]; /* -1.8 */
    k_max = modeldataPtr->parametervector[32]; /* 0.0558 */
    k_min = modeldataPtr->parametervector[33]; /* 0.008 */
    k_abs = modeldataPtr->parametervector[34]; /* 0.057 */
    k_gri = modeldataPtr->parametervector[35]; /* 0.0558 */
    f = modeldataPtr->parametervector[36]; /* 0.9 */
    b = modeldataPtr->parametervector[37]; /* 0.82 */
    d = modeldataPtr->parametervector[38]; /* 0.01 */
    BW = modeldataPtr->parametervector[39]; /* 78 */
    k_p1 = modeldataPtr->parametervector[40]; /* 2.7 */
    k_p2 = modeldataPtr->parametervector[41]; /* 0.0021 */
    k_p3 = modeldataPtr->parametervector[42]; /* 0.009 */
    k_p4 = modeldataPtr->parametervector[43]; /* 0.0618 */
    k_i = modeldataPtr->parametervector[44]; /* 0.0079 */
    U_ii = modeldataPtr->parametervector[45]; /* 1 */
    V_m0 = modeldataPtr->parametervector[46]; /* 2.5 */
    V_mX = modeldataPtr->parametervector[47]; /* 0.047 */
    K_m0 = modeldataPtr->parametervector[48]; /* 225.59 */
    V_f0 = modeldataPtr->parametervector[49]; /* 2.5 */
    V_fX = modeldataPtr->parametervector[50]; /* 0.047 */
    K_f0 = modeldataPtr->parametervector[51]; /* 225.59 */
    p_2U = modeldataPtr->parametervector[52]; /* 0.0331 */
    part = modeldataPtr->parametervector[53]; /* 0.2 */
    K = modeldataPtr->parametervector[54]; /* 2.3 */
    alpha = modeldataPtr->parametervector[55]; /* 0.05 */
    beta = modeldataPtr->parametervector[56]; /* 0.11 */
    gamma = modeldataPtr->parametervector[57]; /* 0.5 */
    k_e1 = modeldataPtr->parametervector[58]; /* 0.0005 */
    k_e2 = modeldataPtr->parametervector[59]; /* 339 */
    D = modeldataPtr->parametervector[60]; /* 78000 */
    /* VARIABLES */
    MeasIR = IR_P;
    MeasIRS = IRS_P;
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
    glucoseuptake = k_glut1*G_t/(KmG1+G_t)+k_glut4Wb*GLUT4_M*G_t/(KmG4+G_t);
    glut1 = k_glut1*G_t/(KmG1+G_t);
    glut4 = k_glut4Wb*GLUT4_M*G_t/(KmG4+G_t);
    v1f = k1f*INS*IR+k1aBasic*IR;
    v1b = k1b*IR_P;
    v2f = k2f*IR_P*IRS;
    v2b = k2b*IRS_P;
    v4f = k4f*IRS_P*PKB;
    v4b = k4b*PKB_P;
    v5f = k5f*PKB_P*GLUT4_C+k5BasicWb*GLUT4_C;
    v5b = k5b*GLUT4_M;
    /* Dependent on the CVODEflag, either the ODERHSvector are determined, */
    /* or the reaction and variable values. */
    if (CVODEflag == CVODE_RHS) {
        /* DIFFERENTIAL EQUATIONS (OUTPUT) */
        ODERHSvector[0] = v1b-v1f;
        ODERHSvector[1] = -v1b+v1f;
        ODERHSvector[2] = v2b-v2f;
        ODERHSvector[3] = -v2b+v2f;
        ODERHSvector[4] = v4b-v4f;
        ODERHSvector[5] = -v4b+v4f;
        ODERHSvector[6] = v5b-v5f;
        ODERHSvector[7] = -v5b+v5f;
        ODERHSvector[8] = EGP+Ra-E-U_ii-k_1*G_p+k_2*G_t;
        ODERHSvector[9] = (-U_id)+k_1*G_p-k_2*G_t;
        ODERHSvector[10] = (-m_1*I_l)-m_3*I_l+m_2*I_p+S;
        ODERHSvector[11] = (-m_2*I_p)-m_4*I_p+m_1*I_l;
        ODERHSvector[12] = -k_gri*Q_sto1;
        ODERHSvector[13] = (-k_empt*Q_sto2)+k_gri*Q_sto1;
        ODERHSvector[14] = (-k_abs*Q_gut)+k_empt*Q_sto2;
        ODERHSvector[15] = -k_i*(I_1-I);
        ODERHSvector[16] = -k_i*(I_d-I_1);
        ODERHSvector[17] = (-p_2U*INS)+p_2U*(I-I_b);
        ODERHSvector[18] = (-gamma*I_po)+S_po;
        ODERHSvector[19] = -alpha*(Y-beta*(G-G_b));
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
        reactionvector[3] = v1f;
        reactionvector[4] = v1b;
        reactionvector[5] = v2f;
        reactionvector[6] = v2b;
        reactionvector[7] = v4f;
        reactionvector[8] = v4b;
        reactionvector[9] = v5f;
        reactionvector[10] = v5b;
    } else if (CVODEflag == CVODE_EVENTS) {
        /* EVENTS */
    } else if (CVODEflag == CVODE_EVENTASSIGNMENT) {
        /* EVENT ASSIGNMENTS */
    }
}
