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
global EXPDATA_IR_TIME
global EXPDATA_IRS1_TIME_30
global EXPDATA_IRS1_TIME_3
global EXPDATA_IRS1_DOUBLE
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
load('EXPDATA_IR_TIME.mat');
load('EXPDATA_IRS1_TIME_30.mat');
load('EXPDATA_IRS1_TIME_3.mat');
load('EXPDATA_IRS1_DOUBLE.mat');
load('EXPDATA_WB.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          THE OPTIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pNamesOpt, startGuess] = SBparameters(optModel);
startGuess = startGuess(3:end);
icOrig = SBinitialconditions(optModel);
icOrigWb = SBinitialconditions(optModelWb);

%startGuess = [11161599367.1635,81973388.5667347,41.2754378139101,0.0187040344795579,0.0343988579029087,3.56363552927496,184.516040553909,1.13952778399062,42.8300150171476,6.32641083472697,0.0954543641709600,696.106720579852,157.119991278241,57.3137982289213,1.27394212321128e-06,6726.02662300000,485.482844100000,0.0209939000000000,0.230788600000000,0.0506481000000000,1.00000000000000e-07,14.7016719000000,0.126915300000000,0.00510630000000000,12485.7473606000,2952.29256510000]';
%   load Result91
%   startGuess = X(1,:)';
startCost = CostFunction(startGuess',0)

OPTIONS = [];
OPTIONS.tempstart = 1e10*startCost;               %InitialTemp
OPTIONS.tempend = 0.5;                           %EndTemp
OPTIONS.tempfactor = 0.1;                      %tempfactor
OPTIONS.maxitertemp = 5000;     %Max-iterations per temp
OPTIONS.maxitertemp0 = 5000;    %Max-iterations at temp0
OPTIONS.maxtime = 2d6;                           %Max-time
OPTIONS.tolx = 1e-10;                            %TolX
OPTIONS.tolfun = 1e-10;                          %Tolfun
OPTIONS.MaxRestartPoints = 3;                    %Number of parallel valleys which are searched through

OPTIONS.lowbounds = [1 1 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9];
OPTIONS.highbounds = [1e12 1e12 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9 1e9];

OPTIONS.outputFunction ='';
OPTIONS.silent = 0;

format long
format compact

% [problem, opts] = ssm_loader(optModel);
% 
% result = ssm_kernel(problem, opts);

[X,FVAL,EXITFLAG] = simannealingSBAOClustering(@CostFunction,startGuess',OPTIONS)


save Result X;

fclose(FID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         PLOT THE RESULT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotData;
plotSimulation(X(1,:));
