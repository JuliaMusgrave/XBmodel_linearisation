%

% This function combines XBmodel_linear with passive_model_linear to provide
% fitting ability of the muslce model
%

% If just the parameters are input then it will output the complex modulus
% using those values,
% if data and freqs are input it will find an objective function such that 
% it can be used in a fitting algorithm to find best parameters
%
% Inputs:
%       - params: vector of fitting paramaters for the model
%           params = [k1 k2 k3i k_1 k_2 K phi_x phi_v phi_s phi_l kPE1 kPE2 eta];
%           no_params = 13
%       - data: 1D array of total complex modulus data to fit to
%       - freqs: 1D array of the sampled frequencies
% Outputs:
%       - OBJ: least-squares error of model vs input data
%       - Y: total complex modulus 
%       - Yp: passive complex modulus
%       - Ya: active complex modulus
%
% Author: Julia Musgrave
% Date: May 2022


function [OBJ,Y,Yp,Ya] = muscle_model_linear(params,data,freqs)

% collecting parameters for each part of the model
xb_params=params(1:10);
p_params=params(11:13);

if nargin ==1
    % solving the parts of the muscle model across default frequency range
    [~,Ya]=XBmodel_linear(xb_params);
    [~,Yp]= passive_model_linear(p_params);
    Y=Ya+Yp;
    OBJ=[];
    return
end

% solving the parts of the muscle model at the input frequencies
[~,Ya]=XBmodel_linear(xb_params,data,freqs);
[~,Yp]=passive_model_linear(p_params,data,freqs);


% Calculating the objective function of the combined model
Y=Ya+Yp;
delY=data-Y;
OBJ=sqrt(0.5/length(data)*sum(delY.*conj(delY)));

% can add a penalty to objective function to force dip
% [dip,idip]=min(abs(Y));
% if (abs(Y(1))-dip)<1e-2 || idip==length(Y)
% OBJ=OBJ+10;
% end

end

