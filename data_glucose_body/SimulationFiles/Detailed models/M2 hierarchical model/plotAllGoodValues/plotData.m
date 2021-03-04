function [] = plotData

clear all
close all

%Figures

figure(10)
h=gca;
hold on
set(h,'box','off');
set(h,'XTick',[0 100 200 300 400])
set(h,'YTick',[0 1 2])
set(h,'Fontsize',[14])
xlabel('time, min')
ylabel('glucose uptake, mg/kg/min')
axis([0 420 0 2])

end

