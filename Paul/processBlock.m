function [I] = processBlock(Y,X,f_s,f_c,dx)
% [I] = processBlock(Y,X,f_s,f_c,dx) take raw data Y and processed for 
% chirp X centered at middle of block. Blocks must be odd in width.
c = 3e8/1.31; %speed of light in ice
n = size(Y,1); %Width
L = size(Y,2); %Time recording length
xx = (((1:n)*dx)-(n+1)/2*dx)'; %lenth vector in azimuth
lambda_c = c/f_c; %Wavelength of center freq of chirp

Imf = ifft( fft(Y,[],2) .* conj(fft(X)) ,[],2);

w = exp(-1i * 2 * pi * (n-1)/2 * [0:floor(n/2)-1 floor(-n/2):-1]'/ n); 
rr = 2 * sqrt(xx.^2 + ((1:L)*(1/f_s)*c/2).^2);
az = exp(1i.*(2.*pi*2*rr./lambda_c));

I = ifft( fft(Imf,[],1)  .* conj( fft(az,[],1) ) .* w,[],1);
