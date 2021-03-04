function [] = plotData

clear all
close all

%Figures
insulinConcentration1=[1e-12,1e-11,1e-10,3e-10,1e-9,1e-8,1e-7]';
insulinConcentration2=[1e-12,1e-11,3e-11,1e-10,3e-10,1e-9,1e-8,1e-7]';

figure(1)
h=gca;
hold on
set(h,'box','off');
set(h,'XTick',[0 100 200 300 400])
set(h,'YTick',[0 100 200])
set(h,'Fontsize',[20])
xlabel('time, min')
ylabel('glucose plasma, mg/dl')
axis([0 420 0 200])

figure(2)
h=gca;
hold on
set(h,'box','off');
set(h,'XTick',[0 100 200 300 400])
set(h,'YTick',[0 100 200 300])
set(h,'Fontsize',[20])
xlabel('time, min')
ylabel('insulin plasma, pM')
axis([0 420 0 300])

figure(10)
h=gca;
hold on
set(h,'box','off');
set(h,'XTick',[0 100 200 300 400])
set(h,'YTick',[0 1 2])
set(h,'Fontsize',[20])
xlabel('time, min')
ylabel('glucose uptake, mg/kg/min')
axis([0 420 0 2])

figure(11)
h=gca;
hold on
set(h,'box','off');
set(h,'XTick',[0 100 200 300 400])
set(h,'YTick',[0 1 2])
set(h,'Fontsize',[20])
xlabel('time, min')
ylabel('IR with multiple insulin bound, %')
axis([0 420 0 2.5])

figure(12)
h=gca;
hold on
set(h,'box','off');
set(h,'XTick',[0 100 200 300 400])
set(h,'YTick',[0 5 10])
set(h,'Fontsize',[20])
xlabel('time, min')
ylabel('PKB phosporylation, %')
axis([0 420 0 11])

end

