clear
%% Inputs 
% Rx
f_s = 10e6;            % Sampling frequency [Hz]
T   = 1e-3;              %Record Time [s];
t = 0:1/f_s:T;       % Time vector [s]
L = length(t);      % Recording length [ ]
% Tx
f_c = 1e6;            % Intial frequency of chirp [Hz] 
swp = 2e6/T;           % Sweep frequency of chirp [Hz/s]
t_c = 1/f_c*100;        % Chirp Length [s]
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

N = length(X);
x = fft(X);
% The mathsy bit. The floors take care of odd-length signals.
w = exp(-1i * 2 * pi * 2*delay * [0:floor(N/2)-1 floor(-N/2):-1] * f_s / N); 
if mod(N, 2) == 0
	% Force conjugate symmetry. Otherwise this frequency component has no
	% corresponding negative frequency to cancel out its imaginary part.
	w(N/2+1) = real(w(N/2+1));
end 
X_r2 = ifft(x .* w * exp(1i.*(4.*pi*r./lambda_c)));

%% Get ready to DFT the pulse, mostly for fun
% This part I took from the internet, I don't really know the details
n = 2^nextpow2(L); % This will pad the signal X with trailing zeros in order to improve the performance of fft.
Y = fft(X,n);

f = f_s*(0:(n/2))/n;
P = abs(Y/n);

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
plot(f,P(1:n/2+1)) 
title('Pulse in Frequency Domain')
xlabel('Frequency (f)')
ylabel('|P(f)|')

disp(mean(angle(MF2).*abs(MF2))/mean(abs(MF2)));

subplot(313)
plot(t,abs((MF1)))
hold on
plot(t,abs((MF2)))
legend('MF1','MF2')
title('Match Filtered Rx in Time Domain')
xlabel('Time (s)')
ylabel('X(t)')