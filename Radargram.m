% This is a script to make n chirps from stationary platform, while moving 
% between chirps to make a scene. Produces both raw data and match filtered
% data products
clear
%% Inputs
% Rx parameters
F_s = 1e7;            % Sampling frequency [Hz]
T   = 3e-4;           % Record Time [s];
t = 0:1/F_s:T;        % Time vector [s]
L = length(t);        % Recording vector length [ ]
% Tx parameters
f_c = 2e6;            % Intial frequency of chirp [Hz] 
swp = 2e6;            % Sweep frequency of chirp [Hz/s]
t_c = 1e-5;           % Chirp Length [s]
% Survey Parameters
n = 101; 			  %surface sample points
dx = 1e3;  			  % Distance between sample points [m]

%% Make pulse
t_sub = 0:1/F_s:t_c; %pulse only for duration of pulse
X = zeros(size(t));  %pulse signal is 0 otherwise
X(1:length(t_sub)) = exp(1i*(pi.*swp.*t_sub.^2+2.*pi.*f_c.*t_sub)); %LFMCM pulse

%% Sweep along surface, let full pulse out and wait for return, then move.
Y = zeros(n,L); %initialize raw data field
Ymf = Y;        %initialize match filtered data field
R = zeros(n,1);
for i = 1:n
    %Find range, shift chirp (in Time and Fx)
    r = sqrt((4e3).^2 + ((i-50)*dx).^2);
    R(i) = r;
    y_tmp = chirpOut(X,t,r,0,f_c);
    
    %raw data
    Y(i,:) = y_tmp; 
    
    % match filter in range, this could be done in a separate processing
    % loop if we wanted, indepentent to data collection
    Ymf(i,:) = ifft( fft(y_tmp)  .* conj( fft(X) ) );
end

%% Plot the returns, abs() of complex values to display
figure(1) %Faster plotting option, but can't see waves
clf
subplot(211)
imagesc(abs(Y'))
ylabel('range')
xlabel('along track')
title('Raw data')
colorbar
subplot(212)
imagesc(abs(Ymf'))
ylabel('range')
xlabel('along track')
title('Match Filtered data')
colorbar

% figure(2) %very slow plotting, but cool to see waveforms!
% clf
% subplot(211)
% wiggle(real(Y'))
% ylabel('range')
% xlabel('along track')
% title('Raw data')
% subplot(212)
% wiggle(abs(Ymf'))
% ylabel('range')
% xlabel('along track')
% title('Range Match Filtered data')
