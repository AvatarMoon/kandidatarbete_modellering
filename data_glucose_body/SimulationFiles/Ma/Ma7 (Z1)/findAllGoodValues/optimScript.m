%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  GLOBALS AND INTIAL STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

global EXPDATA_IR
global EXPDATA_IRS1
global EXPDATA_PKB
global EXPDATA_GLUCOSE_05
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
load('EXPDATA_WB.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          THE OPTIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pNamesOpt, startGuess] = SBparameters(optModel);
icOrig = SBinitialconditions(optModel);
icOrigWb = SBinitialconditions(optModelWb);

startGuess = [1.06856e+006 3344.26 1.38666e+007 659018 1.73005e+006 25432.1 70300.7 614185 864616 1.25831e+007 850673 4.75142e+006 2449.94 498747 477900 40002.9 54329.9 57091.3 77214.2 319.645 365.807 2.44136e+006 68089.9 24244.8 49.4692 128708 12901.9 126969 15863.4 4.31717e+006 783790 1.03402e+007 6267.28 496.969 489.957 599.343 24921.4 0.19725 508689 6.69965e+006 3080.86 1.68558e+006 ]';

startCost = CostFunction(startGuess',0)

OPTIONS = [];
OPTIONS.tempstart = 1e5*startCost;               %InitialTemp
OPTIONS.tempend = 0.5;                           %EndTemp
OPTIONS.tempfactor = 0.1;                        %tempfactor
OPTIONS.maxitertemp = 500*length(startGuess);     %Max-iterations per temp
OPTIONS.maxitertemp0 = 500*length(startGuess);    %Max-iterations at temp0
OPTIONS.maxtime = 2d6;                           %Max-time
OPTIONS.tolx = 1e-10;                            %TolX
OPTIONS.tolfun = 1e-10;                          %Tolfun
OPTIONS.MaxRestartPoints = 10;                    %Number of parallel valleys which are searched through

OPTIONS.lowbounds = [1e-9 1-9 1-9 1-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9];
OPTIONS.highbounds = [1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9];

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
