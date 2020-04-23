function [Y] = chirpOut(X,t,r,theta,f_c,f_s)
% [Y] = chirpOut(X,t,r,theta,f_c) Returns the raw (Y) data
% from an input chirp (X,t,f_c,f_s) off a point scattered range (r) away at angle
% (theta). Chirp (X), time vector(t), center freq (f_c), Sample freq (f_s)
c = 3e8/1.31;

%% Phase delay
lambda_c = c/f_c;
delay = 2*r/c;
%% Attenuation
amp = 1; 
% amp = exp(-r/1e7);

%% Noise, if you want it
SNR = 100;
% noise = randn(size(X)).*exp(1i*rand(size(X))*2*pi); %Rand amp (normal), rand phase (uniform)
noise = zeros(size(X));  %b/c dividing by 0 is hard, use this for 0 noise
%% Shift in space, apply phase and attenuation
%FFT range  shift
N = length(X);
x = fft(X);
% The mathsy bit. The floors take care of odd-length signals.
w = exp(-1i * 2 * pi * delay * [0:floor(N/2)-1 floor(-N/2):-1] * f_s / N); 
if mod(N, 2) == 0
	% Force conjugate symmetry. Otherwise this frequency component has no
	% corresponding negative frequency to cancel out its imaginary part.
	w(N/2+1) = real(w(N/2+1));
end 
Y = ifft(x .* w) * exp(1i*(2*pi*2*r/lambda_c)) + noise./SNR;
% Interpolation range shift
% Y = amp*interp1(t,X,t-delay,'pchip',0).*exp(1i.*(2.*pi*r./lambda_c)) + noise./SNR;