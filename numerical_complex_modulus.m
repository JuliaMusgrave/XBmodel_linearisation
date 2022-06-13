
% This script runs a numerical simulation to determine the complex modulus
% of the ODE-based cross-bridge model XBmodel.m.

% 
% Produces a basic plot of the complex modulus and its components
%
% Author: Julia Musgrave
% Date: May 2022


% input values
model=@XBmodel;
x=[1.825,16.78,6.570,100.0,0.5,9.958e+04,4.177,0.1331,1.133,4.600];
ampl=0.00125; % .125% peak amplitude (.25% pk to pk)
n=25;

% ICs
y0=[0.001 0.001 0 0.0099];

% pre-allocating arrays
Y=zeros(n,1);
Ycomps=zeros(length(y0),n);

% solving ODE system to steady state under constant length
tspan=[0 1];
L0=2.2; % constant length
[~,y]=SSsim(model,tspan,y0,L0,x);
y0=y(end,:);

%% frequency simulations
freqs = logspace(-1,2,n); % n measurements between 0.1 Hz and 100 Hz
for i=1:length(freqs)
f=freqs(i);

T=1/f;
Fs=10000; % sampling frequency
len=Fs*T;

% setting up the t sampling rate and length function
tspan=(0:len-1)/Fs;
L=L0+ampl*L0*sin(2*pi*f*tspan); 

% s is a 2 x len array representing the time-dependent sarcomere length
% input to be applied to the model: row 1 contains time steps, row 2
% contains length values
s=[tspan; L];

% run the sinusoidal sim to steady state
[t,y]=SSsim(model,tspan,y0,s,x);

%calculating force trace
[~,F]=model(t,y,s,x);


% FFT to find complex modulus
YL=fft(L.');
YF=fft(F);

Yf=YF./YL*L0/1000; % Putting into MPa/normalising L
Y(i)=Yf(round(f*len/Fs+1)); %Pulling out complex component at perturbation freq

for j=1:length(y0)
%frequeny response of state variables
Yj=fft(y(:,j));
Yf=Yj./YL*L0/1000;
Ycomps(j,i)=Yf(round(f*len/Fs+1));
end

end

% complex modulus of each state variable scaled as in linearised force eqn
HC_comp = x(6)*y0(4)*Ycomps(2,:);
HxB_comp = x(6)*y0(1)*Ycomps(3,:);
HxC_comp = x(6)*y0(2)*Ycomps(4,:);

%% basic figure of numerical complex modulus
w = 0.4;
h=0.8;

figure('Units', 'normalized' ,'OuterPosition', [0.25, 0.2, 0.45, 0.4])
subplot('Position',[0.075 0.17 w h])
semilogx(freqs,real(Y),'kx','LineWidth',2)
hold on
semilogx(freqs,real(HC_comp),'x','LineWidth',1)
semilogx(freqs,real(HxB_comp),'x','LineWidth',1)
semilogx(freqs,real(HxC_comp),'x','LineWidth',1)
ylabel('Elastic Modulus','FontSize',12)
xlim([0.1 100])
xticklabels({'\fontsize{10} 0.1' '1' '10' '100'})
xlabel('Frequency (Hz)','FontSize',12)
box off 
legend('Full response','HC comp','HxB comp','HxC comp','Location','best')

subplot('Position',[0.575 0.17 w h])
semilogx(freqs,imag(Y),'kx','LineWidth',2)
hold on
semilogx(freqs,imag(HC_comp),'x','LineWidth',1)
semilogx(freqs,imag(HxB_comp),'x','LineWidth',1)
semilogx(freqs,imag(HxC_comp),'x','LineWidth',1)
ylabel('Viscous Modulus','FontSize',12)
xlim([0.1 100])
xticklabels({'\fontsize{10} 0.1' '1' '10' '100'})
xlabel('Frequency (Hz)','FontSize',12)
box off

%%
function [t,y] = SSsim(fun, tspan, y0,s,params)
prev=0;
curr=100;
r=0;
options=odeset('RelTol',1e-6,'Abstol',1e-6,'MaxStep',0.001);

while abs((prev-curr)/curr)>1e-5 && r*tspan(end)<100
[t,y]=ode15s(@(t,y)fun(t,y,s,params),tspan,y0,options);
y0=y(end,:);
prev=curr;
[~,curr]=fun(t(end),y0,s,params);
r=r+1;
end

end