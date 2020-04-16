function [Y] = chirpOut(X,t,r,theta)
% [Y,Y_MF] = chirpOut(X,r) Returns the raw (Y) and Match filtered (Y_MF) data
% from an input chirp (X,t) off a point scattered range (r) away at angle (theta).
delay = r/3e8;

amp = 1; 
% amp = exp(-r/1e7);
Y = amp*interp1(t,X,t-delay,'pchip',0);
