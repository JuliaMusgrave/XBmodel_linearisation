
% This script can be run to plot each of the figures in the paper "Uncovering 
% cross-bridge properties that underlie the cardiac active complex modulus
% using model linearisation techniques"

% Script calls on XBmodel_linear.m, muscle_model_linear.m,
% passive_model_linear.m and param_sensitivity.m functions.
% Also calls on Saeki_data.mat and MLcolours.mat

%
% Author: Julia Musgrave
% Date: May 2022


clear
close all


%% Figure 3 (example of overlap transfer function)
Hov=ones(1,100);
freqs=logspace(-1,2,100);
%[~,~,~,~,Hov_comp] = Land_linear();

plot_plain(freqs,Hov, 'Overlap-based sarcomere length dependence (Fig 3)')

%% Figure 4 (example of state proportion transfer function, with length dependence)
freqs=logspace(-1,2,100);
omi=logspace(-1,2,100)*2*pi*1i;
HC_L=6350./((omi+2095).*(omi+1442)-289.58e+04);

plot_plain(freqs,HC_L, 'State-based sarcomere length dependence (Fig 4)')

%% Figure 5 (example of strain transfer function)
freqs=logspace(-1,2,100);
omi=freqs*2*pi*1i;
HxC=0.83*omi./(omi+29.5);

plot_plain(freqs,HxC,'Sarcomere velocity dependence (Fig 5)')

%% Figure 6 (example of state proportion transfer function, with strain dependence)
freqs=logspace(-1,2,100);
omi=freqs*2*pi*1i;
HxC=0.83*omi./(omi+29.5);
HC_S=(-21*HxC).*(omi+118)./((omi+7).*(omi+118)+22);
plot_plain(freqs,HC_S,'XB strain dependence (Fig 6)')


%% Figure 7 (final fit and component breakdown)

[model,x,YSa,fS]=load_data;
[~,Y,HxB_comp,HxC_comp,HC_comp,HC_Lcomp,HC_Scomp]=model(x);
freqs=logspace(-1,2,100);

[f,A,B,C,D]=four_panel_plot('XB model fit (Fig 7)');

set(f,'CurrentAxes',A)
semilogx(fS,real(YSa),'kx','LineWidth',2,'MarkerSize',10)
semilogx(freqs,real(Y),'k-','LineWidth',2)
text(0.035,max(ylim),'A','FontSize',16,'FontWeight','bold')
l=legend('Active data','Cross-bridge model','Location','northwest','Box','off');
l.ItemTokenSize = [20,18];

set(f,'CurrentAxes',B)
semilogx(fS,imag(YSa),'kx','LineWidth',2,'MarkerSize',10)
semilogx(freqs,imag(Y),'k-','LineWidth',2)
text(0.035,max(ylim),'B','FontSize',16,'FontWeight','bold')

% plotting breakdown of the model
load MLcolours.mat
colors={'k' blue, red, green, green, green};
lines={'-','-','-','-','--','-.'};
comps=[Y;HxB_comp;HxC_comp;HC_comp;HC_Lcomp;HC_Scomp];
width=[2,2,2,2,1.5,1.5];

set(f,'CurrentAxes',C)
for i=1:size(comps,1)
semilogx(freqs,real(comps(i,:)),lines{i},'Color',colors{i},'LineWidth',width(i))
end
text(0.035,max(ylim),'C','FontSize',16,'FontWeight','bold')
l=legend('Full response', '{\it H_{xB}}', '{\it H_{xC}}', '{\it H_C}','length','strain',...
    'Location','northwest','AutoUpdate','off','Box','off');
l.ItemTokenSize = [20,18];
semilogx(freqs,real(Y),'k','LineWidth',2)

set(f,'CurrentAxes',D)
for i=1:size(comps,1)
semilogx(freqs,imag(comps(i,:)),lines{i},'Color',colors{i},'LineWidth',width(i))
end
text(0.035,max(ylim),'D','FontSize',16,'FontWeight','bold')
semilogx(freqs,imag(Y),'k','LineWidth',2)

%% sensitivity analysis plots
[model,x,YSa,fS]=load_data;
param_sensitivity(x,YSa,fS,0);

%Finding the RMSE
RMSE=model(x,YSa,fS);
Ydata=[real(YSa);imag(YSa)];
Range=range(Ydata);
NRMSE=RMSE/Range*100;

%% Figure 9 (comparing breakdown of muscle model fit to total complex modulus)
[model,x,YSa,fS,YSt,YSp]=load_data;
xt=[0.529 105 9.64 91.1	44.2 99700	0.152	0.251	4.88	1.90	183	113	0.326];

[~,Yt,Yp,Ya] = muscle_model_linear(xt);

%finding the RMSE of the XB model when these parameters are used
RMSE_bad=model(xt(1:10),YSa,fS);
Ydata=[real(YSa);imag(YSa)];
Range=range(Ydata);
NRMSE_bad=RMSE_bad/Range*100;

[f,A,B,C,D]=four_panel_plot('Muscle model fit (Fig 9)');

set(f,'CurrentAxes',A)
semilogx(fS,real(YSt),'kx','LineWidth',2,'MarkerSize',8)
semilogx(freqs,real(Yt),'k-','LineWidth',2)
text(0.035,max(ylim),'A','FontSize',16,'FontWeight','bold')
l=legend('Total data','Muscle model','Location','northwest','Box','off');
l.ItemTokenSize = [20,18];

set(f,'CurrentAxes',B)
semilogx(fS,imag(YSt),'kx','LineWidth',2,'MarkerSize',8)
semilogx(freqs,imag(Yt),'k-','LineWidth',2)
text(0.035,max(ylim),'B','FontSize',16,'FontWeight','bold')


% plotting breakdown of model and data into active (XB) and passive
% components
load MLcolours.mat

set(f,'CurrentAxes',C)
semilogx(fS,real(YSa),'x','LineWidth',2,'Color',blue,'MarkerSize',8)
semilogx(freqs,real(Ya),'LineWidth',2,'Color',blue)
semilogx(fS,real(YSp),'x','LineWidth',2,'Color',red,'MarkerSize',8)
semilogx(freqs,real(Yp),'LineWidth',2,'Color',red)
text(0.035,max(ylim),'C','FontSize',16,'FontWeight','bold')
l=legend('Active data', 'Muscle model: Cross-bridge', 'Passive data', 'Muslce model: Passive',...
    'Location','northwest','AutoUpdate','off','Box','off');
l.ItemTokenSize = [20,18];

set(f,'CurrentAxes',D)
semilogx(fS,imag(YSa),'x','LineWidth',2,'Color',blue,'MarkerSize',8)
semilogx(freqs,imag(Ya),'LineWidth',2,'Color',blue)
semilogx(fS,imag(YSp),'x','LineWidth',2,'Color',red,'MarkerSize',8)
semilogx(freqs,imag(Yp),'LineWidth',2,'Color',red)
text(0.035,max(ylim),'D','FontSize',16,'FontWeight','bold')

%% functions

%function to load in Saeki data and model data when needed
function [model,x,Ya,freqs,Yt,Yp]=load_data

load('Saeki_data.mat','Ya','freqs','Yt','Yp')
model=@XBmodel_linear;

%best parameters for model
x=[1.825,16.78,6.570,100.0,0.5,9.958e+04,4.177,0.1331,1.133,4.600];

end

% function to plot EM and VM in a very simple form (for figs 3-6)
function []=plot_plain(fs,Y,name)
w = 0.4;
h=0.72;
ymax=max([real(Y),imag(Y)]);
ymin=min([real(Y),imag(Y)]);
figure('Units', 'normalized' ,'OuterPosition', [0.3, 0.3, 0.27, 0.26],'Name',name)
subplot('Position',[0.075 0.23 w h])
semilogx(fs,real(Y),'k','LineWidth',2)
set(gca,'Fontsize',10)
ylabel('Elastic Modulus','FontSize',12)
xlim([0.1 100])
ylim([min(ymin,0) ymax])
xticklabels({'\fontsize{10} 0.1' '1' '10' '100'})
xlabel('Frequency (Hz)','FontSize',12)
yticks(0)
box off 
text(0.04,max(ylim),'A','FontSize',16,'FontWeight','bold')

subplot('Position',[0.575 0.23 w h])
semilogx(fs,imag(Y),'k','LineWidth',2)
ylabel('Viscous Modulus','FontSize',12)
xlim([0.1 100])
ylim([min(ymin,0) ymax])
xticklabels({'\fontsize{10} 0.1' '1' '10' '100'})
xlabel('Frequency (Hz)','FontSize',12)
yticks(0)
box off
text(0.04,max(ylim),'B','FontSize',16,'FontWeight','bold')

end

% function to set up axes for 4 panel plots (Fig 7 and 9)
function[f,A,B,C,D]=four_panel_plot(name)

f=figure('Name',name, 'Units', 'normalized' ,'OuterPosition', [0.3, 0.3, 0.4, 0.5]);
w = 0.4;
h=0.4;

subplot('Position',[0.075 0.575 w h])
%semilogx(fS,real(YSa),'kx','LineWidth',2,'MarkerSize',10)
hold on
ylabel('Elastic Modulus (MPa)','FontSize',12)
A=gca;
set(A,'Xscale','log','XLim',[0.1 101])
xticks([0.1 1 10 100])
xticklabels({'0.1' '1' '10' '100'})
box off

subplot('Position',[0.575 0.575 w h])
hold on
ylabel('Viscous Modulus (MPa)','FontSize',12)
B=gca;
set(B,'Xscale','log','XLim',[0.1 101])
xticklabels({'0.1' '1' '10' '100'})
box off

subplot('Position',[0.075 0.1 w h])
hold on
ylabel('Elastic Modulus (MPa)','FontSize',12)
C=gca;
set(C,'Xscale','log','XLim',[0.1 101])
xticklabels({'0.1' '1' '10' '100'})
xlabel('Frequency (Hz)','FontSize',12)
box off

subplot('Position',[0.575 0.1 w h])
hold on
ylabel('Viscous Modulus (MPa)','FontSize',12)
D=gca;
set(D,'Xscale','log','XLim',[0.1 101])
xticklabels({'0.1' '1' '10' '100'})
xlabel('Frequency (Hz)','FontSize',12)
box off
end