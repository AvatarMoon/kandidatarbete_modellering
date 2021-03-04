#include "mex.h"

/* Model dependent variables */
const int NRSTATES = 32;
const int NRPARAMETERS = 41;
const int NRVARIABLES = 11;
const int NRREACTIONS = 67;
const int NREVENTS = 0;

double initialconditionsDefault[32] = {
	10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0,10,0};
double parametervectorDefault[41] = {
	0,0,448251462.71204031,4321891.9032703135,0.77226123420000004,0.0122057759,0.0013640432,37.081892484199997,30.682511007700001,0.00095008310000000004,237.51892204340001,3.0181933401999999,0.49752158000000002,128042.884096176,0.0096458629999999993,2374.9773277896002,8.9200000000000005e-008,1.168e-007,608.58395854059995,8.1119350487999995,0.1895302156,180474,1030900,376367,5625.1800000000003,10464600,2380,1859.0699999999999,150488,3348490,35149.699999999997,1.4500000000000001e-008,0.00710824,3.9630700000000001,0.0034996699999999999,1.1900000000000001e-008,59.131799999999998,0.081696900000000003,157.06299999999999,71.620900000000006,0.0037523600000000002};
char *statenames[32] = {
	"r0","r1","r2","r11","r12","r22","r1x2","r11x2","r1x22","r1x22d","r11x22","rend","rendP","iendIR","iend","rPbasal","IRS","IRSiP","X","X_P","PI3K","PI3K_","PDK1","PDK1_","PKC","PKC_P","PKB","PKB_P","mTOR","mTOR_","GLUT4_C","GLUT4_M"};
char *parameternames[41] = {
	"ins","glucose","a1","a2","d1","d2","Kcr","Kex","Kend","Kdp","Kcat","Km","kfbasal","krbasal","k21","k22","k23","k24","k2b","k3f","k3b","k4f","k4b","k5f","k5b","k6f","k6b","k7f","k7b","k8f","k8b","k91","k92","k9b","k5Basic","k5BasicWb","k_glut4","k_glut1","KmG1","KmG4","kbf"};
char *variablenames[11] = {
	"MeasIR","MeasIRS","MeasPKB","MeasAllIntIR","totX","itot","KD","S2","S1","K4","K8"};
char *reactionnames[67] = {
	"vglucoseuptake","vglut1","vglut4","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","R11","R12","R13","R14","R15","R16","R17","R18","R19","R20","R21","R22","R23","R24","R25","R26","R27","R28","R29","R30","R31","R32","R33","R34","R35","R36","R37","R38","R39","R40","R41","R42","R43","R44","R45","R46","R47",
	"R48","v2f","v2b","v3f","v3b","v4f","v4b","v5f","v5b","v6f","v6b","v7f","v7b","v8f","v8b","v9f","v9b"};
char *eventnames[1];

/* Functions containing the model equations */
void model(double time, double *statevector, double *ODERHSvector, ModelData *modeldataPtr, int CVODEflag, double *variablevector, double *reactionvector, double *gout, int *eventvector);

/* CVODE/MEX interface */
void interfaceSER(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    interfaceSER(nlhs, plhs, nrhs, prhs);
}
