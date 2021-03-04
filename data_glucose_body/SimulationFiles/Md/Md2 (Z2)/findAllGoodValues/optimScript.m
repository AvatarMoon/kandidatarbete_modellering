%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  GLOBALS AND INTIAL STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

global EXPDATA_IR
global EXPDATA_IRS1
global EXPDATA_PKB
global EXPDATA_GLUCOSE_05
global EXPDATA_GLUCOSE_5
global EXPDATA_WB
global pNamesOpt  
global icOrig 
global icOrigWb
global modelName
global modelNameWb
global FID

modelName = 'testmodel';
optModel = SBmodel(strcat(modelName,'.txt'));
tmpModel = optModel;
SBAOcreateMEXmodel(optModel,modelName);

modelNameWb = 'testmodelWb';
optModelWb = SBmodel(strcat(modelNameWb,'.txt'));
tmpModelWb = optModelWb;
SBAOcreateMEXmodel(optModelWb,modelNameWb);

FID = fopen('allGoodValues.dat','wt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  CREATE THE EXPDATA STRUCT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('EXPDATA_IR.mat');
load('EXPDATA_IRS1.mat');
load('EXPDATA_PKB.mat');
load('EXPDATA_GLUCOSE_05.mat');
load('EXPDATA_GLUCOSE_5.mat');
load('EXPDATA_WB.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          THE OPTIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pNamesOpt, startGuess] = SBparameters(optModel);
icOrig = SBinitialconditions(optModel);
icOrigWb = SBinitialconditions(optModelWb);

startGuess = [16.3074 20732.7 2.26246e+007 148013 6.41865e+006 274123 361891 389518 56.8715 137559 2.98979e+006 2.22081 2.22081 0.131171 1.55552 1096.59]';
startGuess = 1.0e+008 *[0.000003366995804   0.000238281791382   0.237264545601561  0.001396882175119   0.425470382048549   0.003854564556613   0.000730245229064   0.001235603275933   0.000000088515230   0.000289370836026   0.077022794643419   0.020516190390321   0.098356189843199   0.000000001330966   0.000000001813455   1.607639251170960]';
startGuess = [336.7 23828.2 2.37265e+007 139688 4.2547e+007 385456 73024.5 123560 8.85152 28937.1 7.70228e+006 2.05162e1 9.83562e1 0.133097 0.226682 1.60764e+003 ]';
startGuess = [339.197 24004.9 2.39025e+007 140724 4.28625e+007 388314 72543.1 124476 8.91716 29151.7 7.7594e+006 20.6683 99.0856 0.134084 0.282211 1459.43 ]';

startCost = CostFunction(startGuess',0)

OPTIONS = [];
OPTIONS.tempstart = 1e5*startCost;               %InitialTemp
OPTIONS.tempend = 0.5;                           %EndTemp
OPTIONS.tempfactor = 0.1;                      %tempfactor
OPTIONS.maxitertemp = 500*length(startGuess);     %Max-iterations per temp
OPTIONS.maxitertemp0 = 500*length(startGuess);    %Max-iterations at temp0
OPTIONS.maxtime = 2d6;                           %Max-time
OPTIONS.tolx = 1e-10;                            %TolX
OPTIONS.tolfun = 1e-10;                          %Tolfun
OPTIONS.MaxRestartPoints = 10;                    %Number of parallel valleys which are searched through

OPTIONS.lowbounds = [1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9];
OPTIONS.highbounds = [1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9];

OPTIONS.outputFunction ='';
OPTIONS.silent = 0;

format long
format compact

[X,FVAL,EXITFLAG] = simannealingSBAOClustering(@CostFunction,startGuess',OPTIONS)


save 'Result.m' X;

fclose(FID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         PLOT THE RESULT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotData;
plotSimulation(X(1,:));
