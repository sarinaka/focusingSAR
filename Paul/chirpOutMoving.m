function [Y] = chirpOut(X,t,v,x_0,depth,f_c,f_s)
% [Y] = chirpOut(X,t,v,x_0,depth,theta,f_c,f_s)
% from an input chirp (X,t,f_c,f_s) off a point scattered range (d) and
% inital x offset x_0. Platform moving with velocity v 
% Chirp (X), time vector(t), center freq (f_c), Sample freq (f_s)
c = 3e8/1.79;

r = sqrt((depth).^2 + (x_0 + v*t).^2);
r_c = sqrt((depth).^2 + (x_0).^2);
%% Phase delay
lambda_c = c./f_c;
delay = 2*r./c;
%% Attenuation
amp = 1;
% amp = exp(-r*2/1e3);
%% Noise, if you want it
SNR = 100; %noise ratio to emitted signal
noise = randn(size(X)).*exp(1i*rand(size(X))*2*pi); %Rand amp (normal), rand phase (uniform)
% noise = zeros(size(X));  %b/c dividing by 0 is hard, use this for 0 noise

%% Interp signal forward as option
Rshift = griddedInterpolant(t,X,'spline','none');
% Interpolate signal back up
Y = Rshift(t-delay) .* ... %inter range shift
    exp(1i*(2*pi*(2.*r)/lambda_c)) + ...   % phase delay
    noise./SNR; % noise

Y(isnan(Y))=0;

