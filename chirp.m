clear
%% Inputs 
% Rx
F_s = 1e6;            % Sampling frequency [Hz]
T   = 1;              %Record Time [s];
t = 0:1/F_s:T;       % Time vector [s]
L = length(t);      % Recording length [ ]
% Tx
f_c = 1e3;            % Intial frequency of chirp [Hz] 
swp = 15e5;           % Sweep frequency of chirp [Hz/s]
t_c = 1e-2;            % Chirp Length [s]
% Scatterer
r = 5e6;            %range relay [m]
delay = r/3e8;         %Range delay in time [s]

%% Lets make a chrip!
t_sub = 0:1/F_s:t_c;
X = zeros(size(t));
X(1:length(t_sub)) = exp(1i*(pi.*swp.*t_sub.^2+2.*pi.*f_c.*t_sub));

%% Delay said chirp by a range factor, reflect it back
% Feel free to also reduce by some geometric/attenuation factor
amp = 1;
X_r = chirpOut(X,t,r,0);

%% Get ready to DFT the pulse, mostly for fun
% This part I took from the internet, I don't really know the details
n = 2^nextpow2(L); % This will pad the signal X with trailing zeros in order to improve the performance of fft.
Y = fft(X,n);

f = F_s*(0:(n/2))/n;
P = abs(Y/n);

%% Match Filter 
X_nm = X./(max(abs(X))); %Normalize the ref chirp to 1
MF = ifft( fft(X_r)  .* conj( fft(X) ) ); %should scale ref to 1



figure(1)
clf
subplot(311)
plot(t,real(X))
hold on
plot(t,real(X_r));
title('Pulse in Time Domain')
xlabel('Time (t)')
ylabel('X(t)')
legend('Tx','Rx')

subplot(312)
plot(f,P(1:n/2+1)) 
title('Pulse in Frequency Domain')
xlabel('Frequency (f)')
ylabel('|P(f)|')


subplot(313)
plot(t,abs((MF)))
title('Match Filtered Rx in Time Domain')
xlabel('Time (t)')
ylabel('X(t)')