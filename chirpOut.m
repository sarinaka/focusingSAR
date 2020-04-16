function [Y] = chirpOut(X,t,r,theta,f_c)
% [Y] = chirpOut(X,t,r,theta) Returns the raw (Y) data
% from an input chirp (X,t) off a point scattered range (r) away at angle (theta).
c = 3e8;
lambda_c = c/f_c;
delay = r/c;
amp = 1; 
% amp = exp(-r/1e7);
Y = amp*interp1(t,X,t-delay,'pchip',0).*exp(1i.*(2.*pi*r./lambda_c));
