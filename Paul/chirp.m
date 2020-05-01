clear
%% Inputs 
% Rx
f_s = 10e6;            % Sampling frequency [Hz]
T   = 1e-3;              %Record Time [s];
t = 0:1/f_s:T;       % Time vector [s]
L = length(t);      % Recording length [ ]
% Tx
f_c = 1e5;            % Initial frequency of chirp [Hz] 
swp = 2e6/T;           % Sweep frequency of chirp [Hz/s]
t_c = 1/f_c*10;        % Chirp Length [s]
% Scatterer
c= 3e8/1.31;
r = c./f_c*(100 + pi);   %range relay [m]
delay = r/c;         %Range delay in time [s]
lambda_c = c/f_c;

%% Lets make a chirp!
t_sub = 0:1/f_s:t_c;
X = zeros(size(t));
X(1:length(t_sub)) = exp(1i*(pi.*swp.*t_sub.^2+2.*pi.*f_c.*t_sub));   %Tx

%% Delay said chirp by a range factor, reflect it back
% Feel free to also reduce by some geometric/attenuation factor
amp = 1;
X_r1 = chirpOut(X,t,r,0,f_c,f_s);  %Rx


%% Get ready to DFT the pulse, mostly for fun
% This part I took from the internet, I don't really know the details
n = 2^nextpow2(L); % This will pad the signal X with trailing zeros in order to improve the performance of fft.
Y = fft(X,n);

f = f_s*(0:(n/2))/n;
P = abs(Y/n);

%% Match Filter 
X_nm = X./(max(abs(X))); %Normalize the ref chirp to 1
MF1 = ifft( fft(X_r1)  .* conj( fft(X) ) ); %should scale ref to 1




figure(1)
clf
subplot(311)
plot(t,real(X))
hold on
plot(t,real(X_r1));

title('Pulse in Time Domain')
xlabel('time (s)')
ylabel('X(t)')
legend('Tx','Rx1')

subplot(312)
plot(f,P(1:n/2+1)) 
title('Pulse in Frequency Domain')
xlabel('Frequency (f)')
ylabel('|P(f)|')


subplot(313)
plot(t,abs((MF1)))
legend('MF1')
title('Match Filtered Rx in Time Domain')
xlabel('Time (s)')
ylabel('X(t)')

figure
plot(t,angle(MF1))