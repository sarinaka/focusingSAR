function [I] = processBlock(Y,X,f_s,f_c,dx,t)
% [I] = processBlock(Y,X,f_s,f_c,dx) take raw data Y and processed for 
% chirp X centered at middle of block. Blocks must be odd in width.
c = 3e8/1.31; %speed of light in ice
n = size(Y,1); %Width
L = size(Y,2); %Time recording length
xx = (((1:n)*dx)-(n+1)/2*dx)'; %lenth vector in azimuth
lambda_c = c/f_c; %Wavelength of center freq of chirp

%% Match filter in fast time
Imf = ifft( fft(Y,[],2) .* conj(fft(X)) ,[],2);

%% Need to shift into correct range bin
[x,y] = ndgrid(xx,t*c/2);
Rshift = griddedInterpolant(x,y,Imf);
Imfsh = Rshift(x,t*c/2+x.^2./(2*t*c/2));

%% Match filter in Az
w = exp(-1i * 2 * pi * (n-1)/2 * [0:floor(n/2)-1 floor(-n/2):-1]'/ n); 
rr = sqrt(xx.^2 + ((1:L)*(1/f_s)*c/2).^2);
az = exp(1i.*(2.*pi*(2*rr)./lambda_c));

I = ifft( fft(Imfsh,[],1)  .* conj( fft(az,[],1) ) .* w,[],1);
