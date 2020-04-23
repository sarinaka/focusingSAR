% This is a script to make n chirps from stationary platform, while moving 
% between chirps to make a scene. Produces both raw data and match filtered
% data products
clear
%% Inputs
% Rx parameters
f_s = 40e6;            % Sampling frequency [Hz]
T   = 1e-4;           % Record Time [s];
t = 0:1/f_s:T;        % Time vector [s]
L = length(t);        % Recording vector length [ ]
% Tx parameters
f_c = 5e6;            % Intial frequency of chirp [Hz] 
swp = 15e6/T;            % Sweep frequency of chirp [Hz/s]
t_c = 1e-5;           % Chirp Length [s]
% Survey Parameters
c = 3e8/1.31;
lambda_c = c/f_c;% Wavelength in ice [m]
n = 1001; 			  %surface sample points
dx = 800/1001;  	  % Distance between sample points [m]
depth = 4e3;          % Scatter Depth [m]
% Plotting params
delay_index = round(depth*2/c*f_s);
i_range = 50;
i_min = delay_index-i_range;
i_max = delay_index+i_range;

%% Make pulse
t_sub = 0:1/f_s:t_c; %pulse only for duration of pulse
X = zeros(size(t));  %pulse signal is 0 otherwise
X(1:length(t_sub)) = exp(1i*(pi.*swp.*t_sub.^2+2.*pi.*f_c.*t_sub)); %LFMCM pulse

%% Sweep along surface, let full pulse out and wait for return, then move.
Y = zeros(n,L); %initialize raw data field
Ymf = Y;        %initialize match filtered data field
R = zeros(n,1);
for i = 1:n
    %Find range, shift chirp (in Time and Fx)
    r = sqrt((depth).^2 + ((i-500)*dx).^2);
    R(i) = r;
    y_tmp = chirpOut(X,t,r,0,f_c,f_s);
    
    %raw data
    Y(i,:) = y_tmp; 
    
    % match filter in range, this could be done in a separate processing
    % loop if we wanted, indepentent to data collection
    Ymf(i,:) = ifft( fft(y_tmp)  .* conj( fft(X) ) );
    clear r y_tmp
end

%% Plot the returns, abs() of complex values to display
figure(1) %Faster plotting option, but can't see waves
clf
subplot(311)
	imagesc(real(Y(:,i_min:i_max)'))
	ylabel('range')
	xlabel('along track')
	title('Raw data')
	colorbar
subplot(312)
	imagesc(abs(Ymf(:,i_min:i_max)'))
	ylabel('range')
	xlabel('along track')
	title('Match Filtered data')
	colorbar
subplot(313)
    imagesc(real(Ymf(:,i_min:i_max)'))
    colorbar
    ylabel('range')
    xlabel('along track')
    title('Phase of match filtered data')

% figure(2) %very slow plotting, but cool to see waveforms!
% clf
% subplot(311)
% wiggle(real(Y(1:10:end,:)'))
% ylabel('range')
% xlabel('along track')
% title('Raw data')
% subplot(312)
% wiggle(abs(Ymf(1:10:end,:)'))
% ylabel('range')
% xlabel('along track')
% title('Power of Match Filtered data')
% subplot(313)
% imagesc(angle(Ymf'))
% colorbar
% ylabel('range')
% xlabel('along track')
% title('Phase of Match Filtered data')

%% 
figure(3)
clf
plot(real(Ymf(:,i_min:i_max)))
title('azimuth real data at reflector')

% figure(4)
% clf
% plot((R-depth)/lambda_c)
% title('Range change in units of wavelength')
% 
% figure(5)
% plot(real(exp(1i.*(2.*pi*2*R./lambda_c))))
% title('expected phase shift')