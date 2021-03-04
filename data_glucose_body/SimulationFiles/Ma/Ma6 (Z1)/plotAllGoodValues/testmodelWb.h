#include "mex.h"

/* Model dependent variables */
const int NRSTATES = 25;
const int NRPARAMETERS = 67;
const int NRVARIABLES = 22;
const int NRREACTIONS = 18;
const int NREVENTS = 0;

double initialconditionsDefault[25] = {
	10,0,0,0,0,10,0,10,0,10,0,10,0,178,135,4.5,1.25,78000,0,0,25,25,0,3.6000000000000001,0};
double parametervectorDefault[67] = {
	0,17,0.067511600000000005,0.0092043300000000002,69.209800000000001,0.34086499999999997,0.0036848800000000002,0.118251,0.15049299999999999,1896.1500000000001,4857.4300000000003,858.67100000000005,31102.700000000001,6351.6599999999999,153767,18.196999999999999,1798.8699999999999,928.91700000000003,0.33157500000000001,2.1637,8.9113500000000005,0.112454,1.00417,5100.2200000000003,1.8799999999999999,0.065000000000000002,0.079000000000000001,95,0.050000000000000003,0.19,0.48399999999999999,0.19400000000000001,0.0304,0.64710000000000001,0.59999999999999998,25,1.8,-1.8,0.055800000000000002,0.0080000000000000002,0.057000000000000002,0.055800000000000002,0.90000000000000002,0.81999999999999995,0.01,78,2.7000000000000002,0.0020999999999999999,0.0089999999999999993,0.061800000000000001,
	0.0079000000000000008,1,2.5,0.047,225.59,2.5,0.047,225.59,0.033099999999999997,0.20000000000000001,2.2999999999999998,0.050000000000000003,0.11,0.5,0.00050000000000000001,339,78000};
char *statenames[25] = {
	"IR","IRins","IR_P","IRiP","IRi","IRS","IRSiP","X","X_P","PKB","PKB_P","GLUT4_C","GLUT4_M","G_p","G_t","I_l","I_p","Q_sto1","Q_sto2","Q_gut","I_1","I_d","INS","I_po","Y"};
char *parameternames[67] = {
	"ins","glucose","k1a","k1aBasic","k1b","k1c","k1d","k1e","k1f","k1g","k1r","k21","k22","k2b","k3f","k3b","k4f","k4b","k5f","k5b","k_glut4","k_glut1","KmG1","KmG4","V_G","k_1","k_2","G_b","V_I","m_1","m_2","m_4","m_5","m_6","HE_b","I_b","S_b","S_b_minus","k_max","k_min","k_abs","k_gri","f","b","d","BW","k_p1","k_p2","k_p3","k_p4",
	"k_i","U_ii","V_m0","V_mX","K_m0","V_f0","V_fX","K_f0","p_2U","part","K","alpha","beta","gamma","k_e1","k_e2","D"};
char *variablenames[22] = {
	"MeasIR","MeasIRS","MeasPKB","aa","cc","EGP","V_mmax","V_fmax","E","S","U_idf","I","G","HE","m_3","Q_sto","Ra","k_empt","U_idm","U_id","U","S_po"};
char *reactionnames[18] = {
	"glucoseuptake","glut1","glut4","v1a","v1b","v1c","v1d","v1e","v1g","v1r","v2f","v2b","v3f","v3b","v4f","v4b","v5f","v5b"};
char *eventnames[1];

/* Functions containing the model equations */
void model(double time, double *statevector, double *ODERHSvector, ModelData *modeldataPtr, int CVODEflag, double *variablevector, double *reactionvector, double *gout, int *eventvector);

/* CVODE/MEX interface */
void interfaceSER(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    interfaceSER(nlhs, plhs, nrhs, prhs);
}
