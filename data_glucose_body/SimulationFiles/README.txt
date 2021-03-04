This catalog contains all models tested. It is divided into the three hypotheses Ma, Mb and Mc and the final hierarchical models Mz.
All hypotheses contain several model structures, and the model structures are sometimes tested against more than one dataset.


For every model structure and dataset combination tested you will find a catalog "findAllGoodValues".
The files in that catalog will search the parameter space for acceptable solutions and plot the best solution found.
The catalog "findAllGoodValues" contains:

* the model txt-files "testmodel.txt" and "testmodelWb.txt", where "testmodelWb.txt" contains the whole-body constraints
* an optimization script named "optimScript.m" with changeable optimizing settings
* a function calculating the distance between experimental data and model simulations named "costFunction.m"
* the file "allGoodValues.dat" containing the acceptable parameter sets 
* the scripts "plotData.m" and "plotSimulation.m" for plotting data and simulations
* plots of the best solution (if no acceptable solutions were found), i.e. "IR.eps"
* files containing experimental data, i.e. "EXPDATA_IR.mat"


For the model structure and dataset combinations with acceptable parameter solutions you will also find a catalog "plotAllGoodValues".
The files in that catalog will plot all acceptable parameters. The catalog "plotAllGoodValues" contains:

* the script "plotAllGoodValues.m" that plots all the acceptable parameter sets stored in "allGoodValues.dat"
* the model txt-files "testmodel.txt" and "testmodelWb.txt", where "testmodelWb.txt" contains the whole-body constraints
* the file "allGoodValues.dat" containing the acceptable parameter sets 
* the scripts "plotData.m" and "plotSimulation.m" for plotting data and simulations
* plots of all the acceptable parameters, i.e. "IR.eps"
* files containing experimental data, i.e. "EXPDATA_IR.mat"


To run the m-files SBtoolbox and SBaddon are required (http://www.sbtoolbox.org).
The optimization function "simannealingSBAOClustering.m" has been included in each folder.