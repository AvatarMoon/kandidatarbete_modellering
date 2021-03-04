function[] = plotSimulation(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD THE DATA AND MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('EXPDATA_IR.mat');
load('EXPDATA_IRS1.mat');
load('EXPDATA_PKB.mat');
load('EXPDATA_GLUCOSE_05.mat');
load('EXPDATA_WB.mat');

modelName = 'testmodel';
optModel = SBmodel(strcat(modelName,'.txt'));
tmpModel = optModel;
SBAOcreateMEXmodel(optModel,modelName);

modelNameWb = 'testmodelWb';
optModelWb = SBmodel(strcat(modelNameWb,'.txt'));
tmpModelWb = optModelWb;
SBAOcreateMEXmodel(optModelWb,modelNameWb);

[pNamesOpt, startGuess] = SBparameters(optModel);
icOrig = SBinitialconditions(optModel);
icOrigWb = SBinitialconditions(optModelWb);

names(1:2)=pNamesOpt(1:2)';
names(3:16)=pNamesOpt(3:16)';

simOptions = [];
simOptions.method = 'stiff';
simOptions.maxnumsteps = 1e4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SIMULATE THE MODEL TO STEADY-STATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Simulation to steadystate with mo insulin and no glucose (1) or 130 mg/dl
%glucose (2)
firstparam(1)=0;
firstparam(2)=0;
firstparam(3:16)=param;
firstparamWb(1)=0;
firstparamWb(2)=130;
firstparamWb(3:16)=param;
try
  simDataSteadyState = SBAOsimulate(modelName,1000,icOrig,names,firstparam,simOptions);
  simDataSteadyStateWb = SBAOsimulate(modelName,1000,icOrig,names,firstparamWb,simOptions);
catch disp('Simulation "init" crashed...');
    error = inf;
    return 
end

%New initial conditions
initcond=simDataSteadyState.statevalues(end,:);
initcondWb=simDataSteadyStateWb.statevalues(end,:);
initcondWb(9:20)=icOrigWb(9:20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SIMULATE THE MODEL WITH LOW RESOULTION TO 
%%% CALCULATE THE COST AND TO FIND MAX VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Time-vectors from experimental data
timeGLU=EXPDATA_WB.time;


%Parameter-vectors with different insulin concentrations
for i=1:8
paramP(:,i)=firstparam;
end
paramP(1,:)=[0 0.01 0.1 0.3 1 10 100 0.03];

%Simulations with different insulin concentrations for 15 min (glucose
%uptake) and 10 min (phosphorylations)
try
    for i=1:8
        simData(i) = SBAOsimulate(modelName,15,initcond,names,paramP(:,i),simOptions);
    end
catch disp('Now the P simulation crashed...');
   error = inf;
   return 
end

%Simulation with whole-body insulin and glucose for 420 minutes
try
    simDataGLU = SBAOsimulate(modelNameWb,timeGLU,initcondWb,names,paramP(:,7),simOptions);
catch disp('Now the GLU simulation crashed...');
    error = inf;
    return 
end

%Simulation of dose-response glucose uptake (0.5 and 0.05 mM glucose)
for i=1:8
    initcondGlu(:,i)=simData(i).statevalues(end,:);
    paramGlu05(:,i)=paramP(:,i);
    paramGlu05(2,i)=17;
end
try
    for i=1:8
        simData05(i) = SBAOsimulate(modelName,30,initcondGlu(:,i),names,paramGlu05(:,i),simOptions);
     end
catch disp('Now the Glu dose-response simulation crashed...');
    error = inf;
    return 
end


%Normalization of IR
maxIR = simData(7).variablevalues(667,1);
minIR = simData(1).variablevalues(667,1);
for i=1:7
    simIR(i)= (simData(i).variablevalues(667,1)-minIR)/(maxIR-minIR)*100;
end

%Normalization of IRS1
maxIRS1 = simData(7).variablevalues(667,2);
minIRS1 = simData(1).variablevalues(667,2);
for i=1:7
    simIRS1(i)= (simData(i).variablevalues(667,2)-minIRS1)/(maxIRS1-minIRS1)*100;
end

%Normalization of PKB
maxPKB = simData(7).variablevalues(667,3);
minPKB = simData(1).variablevalues(667,3);
for i=1:7
    simPKB(i)= (simData(i).variablevalues(667,3)-minPKB)/(maxPKB-minPKB)*100;
end

%GLUCOSE UPTAKE 0.5 mM glucose
for i=1:2
    simGLU(i)=simData05(i).reactionvalues(end,1);
    simGLU1(i)=simData05(i).reactionvalues(end,2);
    simGLU4(i)=simData05(i).reactionvalues(end,3);
end
simGLU(3)=simData05(8).reactionvalues(end,1);
simGLU1(3)=simData05(8).reactionvalues(end,2);
simGLU4(3)=simData05(8).reactionvalues(end,3);
for i=3:7
    simGLU(i+1)=simData05(i).reactionvalues(end,1);
    simGLU1(i+1)=simData05(i).reactionvalues(end,2);
    simGLU4(i+1)=simData05(i).reactionvalues(end,3);    
end

%Glucose Uptake WB
for i = 1:11
    simGLUtime(i)=simDataGLU.reactionvalues(i,1);
    simGLUT1(i)=simDataGLU.reactionvalues(i,2);
    simGLUT4(i)=simDataGLU.reactionvalues(i,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATE THE COST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpErrorIR=0;
tmpErrorIRS1=0;
tmpErrorPKB=0;
tmpErrorGLU=0;
tmpErrorWB=0;
for i=1:7
    tmpErrorIR = tmpErrorIR + ((EXPDATA_IR.response(i)-simIR(i))/EXPDATA_IR.sem).^2;
    tmpErrorIRS1 = tmpErrorIRS1 + ((EXPDATA_IRS1.response(i)-simIRS1(i))/EXPDATA_IRS1.sem).^2;
    tmpErrorPKB = tmpErrorPKB + ((EXPDATA_PKB.response(i)-simPKB(i))/EXPDATA_PKB.sem).^2;
end
for i=1:8
    tmpErrorGLU = tmpErrorGLU + ((EXPDATA_GLUCOSE_05.response(i)-simGLU(i))/EXPDATA_GLUCOSE_05.sem).^2;
end
for i=1:11
    tmpErrorWB = tmpErrorWB + ((EXPDATA_WB.response(i)-simGLUtime(i))/EXPDATA_WB.sem).^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SIMULATE WITH HIGH RESOLUTION TO PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timeGLU=0:1:420;

for i=1:8
paramP(:,i)=firstparam;
end
paramP(1,:)=[0 0.01 0.1 0.3 1 10 100 0.03];

%Simulation with whole-body insulin and glucose for 420 minutes
try
    simDataGLU = SBAOsimulate(modelNameWb,timeGLU,initcondWb,names,paramP(:,7),simOptions);
catch disp('Now the simulation Dalla Man crashed...');
    error = inf;
    return 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NORMALIZATION OF DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Glucose Uptake WB
for i = 1:421
    simGLUtime(i)=simDataGLU.reactionvalues(i,1);
    simGLUT1(i)=simDataGLU.reactionvalues(i,2);
    simGLUT4(i)=simDataGLU.reactionvalues(i,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOT THE SIMULATED DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

insulinConcentration1=[1e-12,1e-11,1e-10,3e-10,1e-9,1e-8,1e-7]';
insulinConcentration2=[1e-12,1e-11,3e-11,1e-10,3e-10,1e-9,1e-8,1e-7]';

figure(2)
semilogx(insulinConcentration1,simIR,'b-');
figure(3)
semilogx(insulinConcentration1,simIRS1,'b-');
figure(4)
semilogx(insulinConcentration1,simPKB,'b-');
figure(6)
semilogx(insulinConcentration2,simGLU,'b-');
figure(10)
plot(timeGLU,simGLUtime,'b-');

end