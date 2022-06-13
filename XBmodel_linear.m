
% Function describing a the linearised cross-bridge model presented in "Uncovering 
% cross-bridge properties that underlie the cardiac active complex modulus
% using model linearisation techniques"
%
% Model is fully activated so has no Ca2+ dependence
%
% XB strain equations based on Razumova et al. (1999) 3 state stiffeness 
% distortion model. Sarcomere length dependence on available cross-bridge 
% proportion. Strain dependence on cross-bridge detachment (k3 rate)
% 
% Linearised to calculate the complex modulus under small-amplitude 
% perturbations. 
%
% If just the parameters are input then it will output the complex modulus
% using those values
% if data and freqs are input it will compute the objective function that can be 
% used in a fitting algorithm to find best parameters
%
% Inputs:
%       - params: vector of fitting paramaters for the model
%           params = [k1 k2 k3i k_1 k_2 K phi_x phi_v phi_s phi_l]
%           no_params = 10
%       - data: 1D array of complex modulus to fit to (optional)
%       - freqs: 1D array of the sampled frequencies (optional)
% Outputs:
%       - OBJ: least-squares error of model vs input data (MPa)
%       - Y: active complex modulus (MPa)
%       - HxB_comp: contribution of HxB to the complex modulus (MPa)
%       - HxC_comp: contribution of HxC to the complex modulus (MPa)
%       - HC_comp: contribution of HC to the complex modulus (MPa)
%       - HC_Lcomp: contribution of the length-dependent part of HC to the 
%           complex modulus (MPa)
%       - HC_Scomp: contribution of the strain-dependent part of HC to the 
%           complex modulus (MPa)
%
% Author: Julia Musgrave
% Date: May 2022



function [OBJ,Y,HxB_comp,HxC_comp,HC_comp,HC_Lcomp,HC_Scomp] = XBmodel_linear(params,data,freqs)
if nargin==1
    freqs = logspace(-1,2,100);
end
om=freqs*2*pi;
omi=om*1i; %omega*i

% assigning parameters
k1=params(1); % s^-1
k2=params(2); % s^-1
k3i=params(3); % s^-1
k_1=params(4); % s^-1
k_2=params(5); % s^-1
K=params(6); % kPa
phi_x=params(7); % unitless strain dependence on strain
phi_v=params(8);  % unitless velocity dependence (on strain)
phi_s=params(9); % unitless strain dependence (on k3 rate)
phi_l=params(10); % unitless length dependence (on available XBs)


% model constants
Lmax=2.3; %maximum SL (um)
xC0=0.01;  % powerstroke size (um) 
L0= 2.2; % experiment SL (um)

%----------------------------
% algebraic values at steady state
k3=k3i;  
Z0=1-(Lmax-L0)*phi_l/Lmax;


% thermodynamic constraint for k-3 (as in Tran et al., (2010))
G_ATP0 = -30; % kJ/mol
R = 8.314/1000; % kJ/(mol.K)
T=20+273; % K (Saeki et al. (1991) experiments performed at 20 degC)
ADP=36e-6; Pi=0.001; ATP=0.005; % mol (ATP from Saeki et al. (1991))

G_ATP=G_ATP0+R*T*log((ADP*Pi)/ATP);
k_3 = k1 * k2 * k3 / (k_1 * k_2 * exp(-G_ATP/(R*T)));


% SS state proportions can be found by solving the two ODEs simulataneously
coeffs=[-(k_1+k2+k1), k_2-k1; ...
         k2-k_3, -(k_2+k3+k_3)];
b=[-k1*Z0;-k_3];
SS=coeffs\b;

B0=SS(1);
C0=SS(2);
A0=Z0-B0-C0;

%----------------------------
% partial derivatives
d11=-(k1+k_1+k2);
d21= k_2-k1;
du1= k1*phi_l/Lmax;

d12= k2-k_3;
d22= -(k3+k_2+k_3);
d42= -C0*k3i*phi_s/xC0; 
du2=k_3*phi_l/Lmax;

d33=-phi_x/B0*(A0*k1 + C0*k_2);
du3=phi_v;

d44=-phi_x/C0*(B0*k2+A0*k_3);
du4=phi_v;

%----------------------------
%transfer functions
HxB=du3*omi./(omi-d33);
HxC=du4*omi./(omi-d44);

HC_strain=(d42*HxC.*(omi-d11))./((omi-d22).*(omi-d11)-d12*d21);
HC_len=(du1*d12+du2*(omi-d11))./((omi-d22).*(omi-d11)-d12*d21);
HC=HC_strain+HC_len;

% full response
scale=L0/1000; % to scale CM into MPa (model in kPa)
Y=scale*K*(B0*HxB+C0*HxC+HC*xC0);

% components of full response
HxB_comp=scale*K*B0*HxB;
HxC_comp=scale*K*C0*HxC;
HC_Lcomp=scale*K*xC0*HC_len;
HC_Scomp=scale*K*xC0*HC_strain;
HC_comp=HC_Scomp+HC_Lcomp;


if nargin==1
    OBJ=[];
    return
end

%----------------------------
% calculating RMSE objective function
delY=data-Y;
RMSE=sqrt(0.5/length(data)*sum(delY.*conj(delY)));
OBJ=RMSE; %+(0.032-K*C0*xC0/1000)^2; % option to add steady-state stress to obj

end

