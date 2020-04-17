function [Y] = chirpOut(X,t,r,theta,f_c)
% [Y] = chirpOut(X,t,r,theta,f_c) Returns the raw (Y) data
% from an input chirp (X,t,f_c) off a point scattered range (r) away at angle
% (theta). Chirp (X), time vector(t), center freq (f_c).
c = 3e8/1.31;

%% Phase delay
lambda_c = c/f_c;
delay = r/c;

%% Attenuation
amp = 1; 
% amp = exp(-r/1e7);

%% Noise, if you want it
SNR = 100;
noise = randn(size(X)).*exp(1i*rand(size(X))*2*pi); %Rand amp (normal), rand phase (uniform)
% noise = zeros(size(X));  %b/c dividing by 0 is hard, use this for 0 noise
%% Shift in space, apply phase and attenuation
Y = amp*interp1(t,X,t-delay,'pchip',0).*exp(1i.*(2.*pi*r./lambda_c)) + noise./SNR;
