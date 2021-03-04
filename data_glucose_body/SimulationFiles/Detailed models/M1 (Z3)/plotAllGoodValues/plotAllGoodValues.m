function [] = plotAllGoodValues()

A=load('allGoodValues.dat');

[sizeGoodParamSets sizeParam]=size(A);

if sizeGoodParamSets>0
    % Find max and min values of all parameters
    extremevalues = [];
    for i=1:sizeParam
        extremevalues(i,1)=max(A(:,i));
        extremevalues(i,2)=min(A(:,i));
    end
    
    % Find the positions in the matrix where max and min is situated
    for i=1:sizeParam
        position(i,1)=find(A(:,i)>=extremevalues(i,1),1);
        position(i,2)=find(A(:,i)<=extremevalues(i,2),1);
    end
    % Plot experimental data
    plotData

    %Plot simulated data (min and max parameters)
    for i=1:sizeParam
        for j=1:2
            param=A(position(i,j),2:end);
            plotSimulation(param)
        end
    end
else
    disp('There are no acceptable solutions!')
end
end