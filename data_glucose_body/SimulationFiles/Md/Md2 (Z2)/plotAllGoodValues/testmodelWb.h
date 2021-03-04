#include "mex.h"

/* Model dependent variables */
const int NRSTATES = 20;
const int NRPARAMETERS = 61;
const int NRVARIABLES = 22;
const int NRREACTIONS = 11;
const int NREVENTS = 0;

double initialconditionsDefault[20] = {
	10,0,10,0,10,0,10,0,178,135,4.5,1.25,78000,0,0,25,25,0,3.6000000000000001,0};
double parametervectorDefault[61] = {
	0,17,107610,720.74599999999998,726898,479.07799999999997,999792,119471,159236,0.001,0.001,30433.900000000001,47736.800000000003,79749.300000000003,79749.300000000003,1.6000000000000001,294,21,1.8799999999999999,0.065000000000000002,0.079000000000000001,95,0.050000000000000003,0.19,0.48399999999999999,0.19400000000000001,0.0304,0.64710000000000001,0.59999999999999998,25,1.8,-1.8,0.055800000000000002,0.0080000000000000002,0.057000000000000002,0.055800000000000002,0.90000000000000002,0.81999999999999995,0.01,78,2.7000000000000002,0.0020999999999999999,0.0089999999999999993,0.061800000000000001,0.0079000000000000008,1,2.5,0.047,225.59,2.5,
	0.047,225.59,0.033099999999999997,0.20000000000000001,2.2999999999999998,0.050000000000000003,0.11,0.5,0.00050000000000000001,339,78000};
char *statenames[20] = {
	"IR","IR_P","IRS","IRS_P","PKB","PKB_P","GLUT4_C","GLUT4_M","G_p","G_t","I_l","I_p","Q_sto1","Q_sto2","Q_gut","I_1","I_d","INS","I_po","Y"};
char *parameternames[61] = {
	"ins","glucose","k1aBasic","k1f","k1b","k2f","k2b","k4f","k4b","k5Basic","k5BasicWb","k5f","k5b","k_glut4","k_glut4Wb","k_glut1","KmG1","KmG4","V_G","k_1","k_2","G_b","V_I","m_1","m_2","m_4","m_5","m_6","HE_b","I_b","S_b","S_b_minus","k_max","k_min","k_abs","k_gri","f","b","d","BW","k_p1","k_p2","k_p3","k_p4","k_i","U_ii","V_m0","V_mX","K_m0","V_f0",
	"V_fX","K_f0","p_2U","part","K","alpha","beta","gamma","k_e1","k_e2","D"};
char *variablenames[22] = {
	"MeasIR","MeasIRS","MeasPKB","aa","cc","EGP","V_mmax","V_fmax","E","S","U_idf","I","G","HE","m_3","Q_sto","Ra","k_empt","U_idm","U_id","U","S_po"};
char *reactionnames[11] = {
	"glucoseuptake","glut1","glut4","v1f","v1b","v2f","v2b","v4f","v4b","v5f","v5b"};
char *eventnames[1];

/* Functions containing the model equations */
void model(double time, double *statevector, double *ODERHSvector, ModelData *modeldataPtr, int CVODEflag, double *variablevector, double *reactionvector, double *gout, int *eventvector);

/* CVODE/MEX interface */
void interfaceSER(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    interfaceSER(nlhs, plhs, nrhs, prhs);
}
