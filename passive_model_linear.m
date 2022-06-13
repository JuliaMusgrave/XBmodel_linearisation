%

% A linearised model of passive force based on Tewari et al. (2016) 
%
% If just the parameters are input then it will output the complex modulus
% using those values,
% if data and freqs and are input it will find the complex modulus at the
% specified freqs and calculate the RMSE
%
% Inputs:
%       - params: vector of fitting paramaters for the model
%           params = [kPE1 kPE2 eta]
%           no_params = 3
%       - data: 1D array of complex modulus to fit to
%       - freqs: 1D array of the sampled frequencies
% Outputs:
%       - OBJ: least-squares error of model vs input data
%       - Y: passive complex modulus 
%
% Author: Julia Musgrave
% Date: May 2022


function [OBJ,Y] = passive_model_linear(params,data,freqs)
if nargin==1
    freqs = logspace(-1,2,100);
end

% model parameters
kPE1=params(1);
kPE2=params(2);
eta=params(3);

% frequency response parameters
om=freqs*2*pi;
omi=om*1i; %omega*i
L0= 2.2; % experiment SL (um)

% partial derivatives
d11=-kPE1/eta;
du1=kPE1;

HF1=du1*omi./(omi-d11);

% passive complex modulus
Y=L0/1000*(HF1+kPE2);

if nargin==1
    OBJ=[];
    return
end

delY=data-Y;
OBJ=sqrt(0.5/length(data)*sum(delY.*conj(delY)));


end

