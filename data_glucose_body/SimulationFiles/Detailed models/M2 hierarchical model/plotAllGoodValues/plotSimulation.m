function[] = plotSimulation(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOAD THE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelName = 'hierarchicalmodel';
optModel = SBmodel(strcat(modelName,'.txt'));
tmpModel = optModel;
SBAOcreateMEXmodel(optModel,modelName);

[pNamesOpt, startGuess] = SBparameters(optModel);
icOrig = SBinitialconditions(optModel);

names=pNamesOpt(1:36)';

simOptions = [];
simOptions.method = 'stiff';
simOptions.maxnumsteps = 1e4;

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

%Glucose Uptake
for i = 1:421
    simGluUpt(i)=simData.variablevalues(i,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOT THE SIMULATED DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10)
plot(time,simGluUpt,'b-');

end