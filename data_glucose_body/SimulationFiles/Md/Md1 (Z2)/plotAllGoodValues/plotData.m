function [] = plotData

clear all
close all

%load data
load('EXPDATA_IR.mat');
load('EXPDATA_IRS1.mat');
load('EXPDATA_PKB.mat');
load('EXPDATA_GLUCOSE_05.mat');
load('EXPDATA_GLUCOSE_5.mat');
load('EXPDATA_WB.mat');

%Figures
insulinConcentration1=[1e-12,1e-11,1e-10,3e-10,1e-9,1e-8,1e-7]';
insulinConcentration2=[1e-12,1e-11,3e-11,1e-10,3e-10,1e-9,1e-8,1e-7]';

figure(2)
h=gca;
hold on
for i=1:7
    plot(insulinConcentration1(i),EXPDATA_IR.response(i),'ro','MarkerFaceColor','r','MarkerSize',10,'Linewidth',2)
    line([insulinConcentration1(i) insulinConcentration1(i)],[EXPDATA_IR.response(i)-EXPDATA_IR.sem EXPDATA_IR.response(i)+EXPDATA_IR.sem],'Color','r','Linewidth',2)
end
set(h,'box','off')
set(h,'XScale','log')
set(h,'XTick',[1e-010 1e-008 1e-006])
set(h,'YTick',[0 50 100])
set(h,'Fontsize',[14])
xlabel('[insulin], mol/L')
ylabel('IR phosphorylation, % of max')
htext=text(9e-13,-16,'0');
set(htext,'Fontsize',[14])
axis([7e-013 1.5e-007 -10 110])

figure(3)
h=gca;
hold on
for i=1:7
    plot(insulinConcentration1(i),EXPDATA_IRS1.response(i),'ro','MarkerFaceColor','r','MarkerSize',10,'Linewidth',2)
    line([insulinConcentration1(i) insulinConcentration1(i)],[EXPDATA_IRS1.response(i)-EXPDATA_IRS1.sem EXPDATA_IRS1.response(i)+EXPDATA_IRS1.sem],'Color','r','Linewidth',2)
end
set(h,'box','off')
set(h,'XScale','log')
set(h,'XTick',[1e-010 1e-008 1e-006])
set(h,'YTick',[0 50 100])
set(h,'Fontsize',[14])
xlabel('[insulin], mol/L')
ylabel('IRS1 phosphorylation, % of max')
htext=text(9e-13,-16,'0');
set(htext,'Fontsize',[14])
axis([7e-013 1.5e-007 -10 110])

figure(4)
h=gca;
hold on
for i=1:7
    plot(insulinConcentration1(i),EXPDATA_PKB.response(i),'ro','MarkerFaceColor','r','MarkerSize',10,'Linewidth',2)
    line([insulinConcentration1(i) insulinConcentration1(i)],[EXPDATA_PKB.response(i)-EXPDATA_PKB.sem EXPDATA_PKB.response(i)+EXPDATA_PKB.sem],'Color','r','Linewidth',2)
end
set(h,'box','off')
set(h,'XScale','log')
set(h,'XTick',[1e-010 1e-008 1e-006])
set(h,'YTick',[0 50 100])
set(h,'Fontsize',[14])
xlabel('[insulin], mol/L')
ylabel('PKB phosphorylation, % of max')
htext=text(9e-13,-16,'0');
set(htext,'Fontsize',[14])
axis([7e-013 1.5e-007 -10 110])

figure(6)
h=gca;
hold on
for i=1:8
    plot(insulinConcentration2(i),EXPDATA_GLUCOSE_05.response(i),'ro','MarkerFaceColor','r','MarkerSize',10,'Linewidth',2)
    line([insulinConcentration2(i) insulinConcentration2(i)],[EXPDATA_GLUCOSE_05.response(i)-EXPDATA_GLUCOSE_05.sem EXPDATA_GLUCOSE_05.response(i)+EXPDATA_GLUCOSE_05.sem],'Color','r','Linewidth',2)
end
set(h,'box','off')
set(h,'XScale','log')
set(h,'XTick',[1e-011 1e-009 1e-007])
set(h,'YTick',[0 0.1 0.2 0.3])
set(h,'Fontsize',[14])
xlabel('[insulin], mol/L')
ylabel('glucose uptake, mg/kg/min')
htext=text(9e-13,-0.015,'0');
set(htext,'Fontsize',[14])
axis([7e-013 1.5e-007 0 0.28])

figure(7)
h=gca;
hold on
for i=1:2
    plot(insulinConcentration2(7*i-6),EXPDATA_GLUCOSE_5.response(i),'ro','MarkerFaceColor','r','MarkerSize',10,'Linewidth',2)
    line([insulinConcentration2(7*i-6) insulinConcentration2(7*i-6)],[EXPDATA_GLUCOSE_5.response(i)-EXPDATA_GLUCOSE_5.sem(i) EXPDATA_GLUCOSE_5.response(i)+EXPDATA_GLUCOSE_5.sem(i)],'Color','r','Linewidth',2)
end
set(h,'box','off')
set(h,'XScale','log')
set(h,'XTick',[1e-011 1e-009 1e-007])
set(h,'YTick',[0 1 2])
set(h,'Fontsize',[14])
xlabel('[insulin], mol/L')
ylabel('glucose uptake, mg/kg/min')
htext=text(9e-13,-0.09,'0');
set(htext,'Fontsize',[14])
axis([7e-013 1.5e-007 0 1.5])

figure(10)
h=gca;
hold on
for i=1:11
    plot(EXPDATA_WB.time(i),EXPDATA_WB.response(i),'ro','MarkerFaceColor','r','MarkerSize',10,'Linewidth',2);
    line([EXPDATA_WB.time(i) EXPDATA_WB.time(i)],[EXPDATA_WB.response(i)-EXPDATA_WB.sem EXPDATA_WB.response(i)+EXPDATA_WB.sem],'Color','r','Linewidth',2)
end
set(h,'box','off');
set(h,'XTick',[0 100 200 300 400])
set(h,'YTick',[0 0.5 1])
set(h,'Fontsize',[14])
xlabel('time, min')
ylabel('glucose uptake, mg/kg/min')
axis([0 420 0 1.1])

end
