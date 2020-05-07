clear
%% Inputs 
% Rx
f_s = 10e6;            % Sampling frequency [Hz]
T   = 2e-4;              %Record Time [s];
t = 0:1/f_s:T;       % Time vector [s]
L = length(t);      % Recording length [ ]

% Tx
f_c = 1e5;            % Initial frequency of chirp [Hz] 
t_c = 1/f_c*10;        % Chirp Length [s]
swp = 1e6/t_c;           % Sweep frequency of chirp [Hz/s]
v = 1000000;               %platform velocity [m/s]

% Scatterer
c = 3e8/1.31;
d = 0e3;
x_0 = 2e3;
r1 = sqrt((x_0).^2 + d^2);   %range relay [m]
r2 = sqrt((x_0-v*t).^2 + d^2);   %range relay [m]
disp( "movement in chirp:" + (r2(end)-r2(1)))
lambda_c = c/f_c;

%% Lets make a chirp!
t_sub = 0:1/f_s:t_c;
X = zeros(size(t));
X(1:length(t_sub)) = exp(1i*(pi.*swp.*t_sub.^2+2.*pi.*f_c.*t_sub));   %Tx

%% Delay said chirp by a range factor, reflect it back
% Feel free to also reduce by some geometric/attenuation factor
amp = 1;
X_r1 = chirpOut(X,t,r1,0,f_c,f_s);  %Rx
X_r2 = chirpOutMoving(X,t,v,x_0,d,f_c,f_s);  %Rx

%% Get ready to DFT the pulse, mostly for fun
% This part I took from the internet, I don't really know the details
n = 2^nextpow2(L); % This will pad the signal X with trailing zeros in order to improve the performance of fft.
Y1 = fft(X_r1,n);
Y2 = fft(X_r2,n);

f = f_s*(0:(n/2))/n;
P1 = abs(Y1/n);
P2 = abs(Y2/n);

%% Match Filter 
X_nm = X./(max(abs(X))); %Normalize the ref chirp to 1
MF1 = ifft( fft(X_r1)  .* conj( fft(X) ) ); %should scale ref to 1
MF2 = ifft( fft(X_r2)  .* conj( fft(X) ) ); %should scale ref to 1

figure(1)
clf
subplot(311)
plot(t,real(X))
hold on
plot(t,real(X_r1));
plot(t,real(X_r2));

title('Pulse in Time Domain')
xlabel('time (s)')
ylabel('X(t)')
legend('Tx','Rx1','Rx2')

subplot(312)
plot(f,P1(1:n/2+1))
hold on
plot(f,P2(1:n/2+1)) 
title('Pulse in Frequency Domain')
xlabel('Frequency (f)')
ylabel('|P|')
legend('|P(Rx1)|','|P(Rx2)|')

subplot(313)
plot(t,abs((MF1)))
hold on
plot(t,abs((MF2)))
legend('MF1','MF2')
title('Match Filtered Rx in Time Domain')
xlabel('Time (s)')
ylabel('X(t)')
