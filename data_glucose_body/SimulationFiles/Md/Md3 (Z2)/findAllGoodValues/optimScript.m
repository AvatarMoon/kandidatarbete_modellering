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

startGuess = [1863.78 38668.3 4.0229e+007 401394 3.67043e+007 336782 187399 85530.9 6.01237 11264.8 2.63899e+007 0.00394328 19.7881 0.255699 172.027 74.5341 ]';

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

OPTIONS.lowbounds = [1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1 1];
OPTIONS.highbounds = [1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e4 1e4];

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
