function[] = plotSimulation(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD THE MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelName = 'hierarchicalmodel';
optModel = SBmodel(strcat(modelName,'.txt'));
tmpModel = optModel;
SBAOcreateMEXmodel(optModel,modelName);

[pNamesOpt, startGuess] = SBparameters(optModel);
icOrig = SBinitialconditions(optModel);

names = pNamesOpt(1:39);

simOptions = [];
simOptions.method = 'stiff';
simOptions.maxnumsteps = 1e6;
simOptions.reltol = 1e-6;
simOptions.abstol = 1e-9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SIMULATE WITH HIGH RESOLUTION TO PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time=0:1:420;

%Simulation with whole-body insulin and glucose for 420 minutes
try
    simData = SBAOsimulate(modelName,time,icOrig,names,param,simOptions);
catch disp('Now the simulation Dalla Man crashed...');
    error = inf;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NORMALIZATION OF DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Glucose plasma
for i = 1:421
    simGluPla(i)=simData.variablevalues(i,14);
end

%Insulin plasma
for i = 1:421
    simInsPla(i)=simData.variablevalues(i,13);
end

%Glucose Uptake
for i = 1:421
    simGluUpt(i)=simData.reactionvalues(i,1);
end

%IR with more than one insulin molecule bound
for i = 1:421
    simIR(i)=0;
    for j = 4:6
        simIR(i)=simIR(i)+(simData.statevalues(i,j))*10;
    end
    for j = 8:11
        simIR(i)=simIR(i)+(simData.statevalues(i,j))*10;
    end
end

%PKB state
for i = 1:421
    simPKB(i)=simData.statevalues(i,28)*10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOT THE SIMULATED DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(time,simGluPla,'b-');
figure(2)
plot(time,simInsPla,'b-');
figure(10)
plot(time,simGluUpt,'b-');
figure(11)
plot(time,simIR,'b-');
figure(12)
plot(time,simPKB,'b-');

end