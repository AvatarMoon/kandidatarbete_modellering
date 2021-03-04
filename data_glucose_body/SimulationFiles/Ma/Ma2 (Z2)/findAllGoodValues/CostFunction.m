function [error] = CostFunction(param,shouldIPlot)

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

simOptions = [];
simOptions.method = 'stiff';
simOptions.maxnumsteps = 1e4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SIMULATE THE MODEL FOR THE GIVEN param
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Selection of parameters to optimize
names(1:2)=pNamesOpt(1:2)';
names(3:16)=pNamesOpt(3:16)';

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

%Time-vectors from experimental data
timeGLU=EXPDATA_WB.time;

%Parameter-vectors with different insulin concentrations
for i=1:8
paramP(:,i)=firstparam;
end
paramP(1,:)=[0 0.01 0.1 0.3 1 10 100 0.03];

%Simulations with different insulin concentrations for 15 min, 20 min
%(glucose uptake) and 10 min (phosphorylations)
try
    for i=1:8
        simData(i) = SBAOsimulate(modelName,20,initcond,names,paramP(:,i),simOptions);
    end
catch disp('Now the P simulation crashed...');
   error = inf;
   return 
end

%Simulation with whole-body insulin and glucose for 420 minutes
try
    simDataGLU = SBAOsimulate(modelNameWb,timeGLU,initcondWb,names,paramP(:,1),simOptions);
catch disp('Now the GLU simulation crashed...');
    error = inf;
    return 
end

%Simulation of dose-response glucose uptake (0.5 mM glucose)
for i=1:8
    initcondGlu(:,i)=simData(i).statevalues(750,:);
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

%Simulation of dose-response glucose uptake (5 mM glucose)
initcondGlu5(:,1)=simData(1).statevalues(1000,:);
initcondGlu5(:,2)=simData(7).statevalues(1000,:);
paramGlu5(:,1)=paramP(:,1);
paramGlu5(:,2)=paramP(:,7);
paramGlu5(2,1)=170;
paramGlu5(2,2)=170;
try
    for i=1:2
        simData5(i) = SBAOsimulate(modelName,8,initcondGlu5(:,i),names,paramGlu5(:,i),simOptions);
    end
catch disp('Now the Glu dose-response simulation crashed...');
    error = inf;
    return 
end

%Normalization of IR
maxIR = simData(7).variablevalues(500,1);
minIR = simData(1).variablevalues(500,1);
for i=1:7
    if maxIR > minIR
        simIR(i)= (simData(i).variablevalues(500,1)-minIR)/(maxIR-minIR)*100;
    else
        simIR(i) = -1000;
    end
end

%Normalization of IRS1
maxIRS1 = simData(7).variablevalues(500,2);
minIRS1 = simData(1).variablevalues(500,2);
for i=1:7
    if maxIRS1 > minIRS1
        simIRS1(i)= (simData(i).variablevalues(500,2)-minIRS1)/(maxIRS1-minIRS1)*100;
    else
        simIRS1(i) = -1000;
    end
end

%Normalization of PKB
maxPKB = simData(7).variablevalues(500,3);
minPKB = simData(1).variablevalues(500,3);
for i=1:7
    if maxPKB > minPKB
        simPKB(i)= (simData(i).variablevalues(500,3)-minPKB)/(maxPKB-minPKB)*100;
    else
        simPKB(i) = -1000;
    end
end

%GLUCOSE UPTAKE 0.5 mM glucose
for i=1:2
    simGLU(i)=simData05(i).reactionvalues(end,1);
end
simGLU(3)=simData05(8).reactionvalues(end,1);
for i=3:7
    simGLU(i+1)=simData05(i).reactionvalues(end,1);    
end

%GLUCOSE UPTAKE 5.0 mM glucose
for i=1:2
    simGLU5(i)=simData5(i).reactionvalues(end,1);
end

%Glucose Uptake WB
for i = 1:11
    simGLUtime(i)=simDataGLU.reactionvalues(i,1);
    simGLUT1(i)=simDataGLU.reactionvalues(i,2);
    simGLUT4(i)=simDataGLU.reactionvalues(i,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       CALCULATE THE COST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmpError = 0;
for i=1:7
   tmpError = tmpError + ((EXPDATA_IR.response(i)-simIR(i))/EXPDATA_IR.sem).^2;
   tmpError = tmpError + ((EXPDATA_IRS1.response(i)-simIRS1(i))/EXPDATA_IRS1.sem).^2;
   tmpError = tmpError + ((EXPDATA_PKB.response(i)-simPKB(i))/EXPDATA_PKB.sem).^2;
end
for i=1:8
    tmpError = tmpError + ((EXPDATA_GLUCOSE_05.response(i)-simGLU(i))/EXPDATA_GLUCOSE_05.sem).^2;
end
for i=1:2
    tmpError = tmpError + ((EXPDATA_GLUCOSE_5.response(i)-simGLU5(i))/EXPDATA_GLUCOSE_5.sem(i)).^2;
end
for i=1:11
    tmpError = tmpError + ((EXPDATA_WB.response(i)-simGLUtime(i))/EXPDATA_WB.sem).^2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    FINAL CHECKS AND RETURN COMMAND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if tmpError < 50
    fprintf(FID, '%4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f \n',[tmpError, param]);
  end
%And this is the return command
  error = tmpError; 
end
