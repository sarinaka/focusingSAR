% This is a script to make n chirps from stationary platform, while moving 
% between chirps to make a scene. Produces both raw data and match filtered
% data products
clear
%% Inputs
% Rx parameters
f_s = 60e6;          % Sampling frequency [Hz]
T   = 5e-5;           % Record Time [s];
t = 0:1/f_s:T;        % Time vector [s]
L = length(t);        % Recording vector length [ ]
% Tx parameters
f_0 = 20e6;           % Intial frequency of chirp [Hz] 
swp = 40e6;           % BW of chirp [Hz]
f_c = f_0 + swp/2;    % Center frequency [Hz]
           
t_c = 1e-5;           % Chirp Length [s]
% Survey Parameters
c = 3e8/1.31;
lambda_c = c/f_c;     % Wavelength (center) in ice [m]
lambda_0 = c/f_0;     % Wavelength (initial) in ice [m]
n = 1001; 			  % surface sample points
dx = 800/1001;  	  % Distance between sample points [m]
dy = (1/f_s)*c/2;     % Distance of 1 way dist per range bin [m]
SynthAx = n*dx;       % Synthetic Apeture [m]
xx = (([1:n]*dx)-(n+1)/2*dx)';
depth = 4e3;          % Scatter Depth [m]
% Plotting params
delay_index = round(depth*2/c*f_s);
i_range = 100;
i_min = delay_index-i_range;
i_max = delay_index+i_range;

%% Make pulse
t_sub = 0:1/f_s:t_c; %pulse only for duration of pulse
X = zeros(size(t));  %pulse signal is 0 otherwise
X(1:length(t_sub)) = exp(1i*(pi.*(swp/t_c).*t_sub.^2+2.*pi.*f_c.*t_sub)); %LFMCM pulse

%% Sweep along surface, let full pulse out and wait for return, then move.
Y = zeros(n,L); %initialize raw data field
Ymf = Y;        %initialize match filtered data field
R = zeros(n,1);
for i = 1:n
    %Find range, shift chirp (in Time and Fx)
    r1 = sqrt((depth).^2 + ((i-500)*dx).^2);
    r2 = sqrt((depth+4e1).^2 + ((i-450)*dx).^2);
    y_tmp = chirpOut(X,t,r1,0,f_c,f_s) + chirpOut(X,t,r2,0,f_c,f_s);
    
    %raw data
    Y(i,:) = y_tmp; 
    
    % match filter in range, this could be done in a separate processing
    % loop if we wanted, indepentent to data collection
    Ymf(i,:) = ifft( fft(y_tmp)  .* conj( fft(X) ) );
    clear r y_tmp
end
clear i
tic
for i = 1:n
    % match filter in range, this could be done in a separate processing
    % loop if we wanted, indepentent to data collection
    Ymf(i,:) = ifft( fft(Y(i,:))  .* conj( fft(X) ) );
    clear r y_tmp
end

[x,y] = ndgrid(xx,t*c/2);
Rshift = griddedInterpolant(x,y,Ymf);
Ymfsh = Rshift(x,t*c/2+x.^2./(2*t*c/2));



%% Azimuth Processing
Ymfaz = zeros(n,i_range*2+1);
%shift in reflector to our domain
w = exp(-1i * 2 * pi * (n-1)/2 * [0:floor(n/2)-1 floor(-n/2):-1]'/ n); 
for j = 1:L
    rr = sqrt(xx.^2 + ((j)*(1/f_s)*c/2)^2);
    z = Ymfsh(:,j);
    az = exp(1i.*(2.*pi*(2*rr)./lambda_c));
    Ymfaz(:,j) = ifft( fft(z)  .* conj( fft(az) ) .* w);
end

toc

I = processBlock(Y,X,f_s,f_c,dx,t); %Processing Az with initial freq (not center) works way better

%% Plot the returns, abs() of complex values to display
figure(1) %Faster plotting option, but can't see waves
clf
subplot(311)
	prettyPlot(real(Y(:,i_min:i_max)'))
	ylabel('range')
	xlabel('along track')
	title('Raw data')
	colorbar
subplot(312)
	prettyPlot(abs(Ymfsh(:,i_min:i_max)'))
	ylabel('range')
	xlabel('along track')
	title('Match Filtered data')
	colorbar
subplot(313)
    prettyPlot(abs(I(:,i_min:i_max)'))
%     surf(abs(Ymfaz(:,i_min:i_max)'),'edgecolor','none')
    colorbar
    ylabel('range')
    xlabel('along track')
    title('Azimuth Focused Data')

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

%% Azimuth Focusing play area

i_max = delay_index + 3;
z = Ymf(:,delay_index:i_max);

% Generate 2-way distance to reflector at this bin
rr = sqrt(xx.^2 + ((delay_index:i_max)*(1/f_s)*c/2).^2);
% Generate phase delay at this distance along Az to filter with
az = exp(1i.*(2.*pi*(2*rr)./lambda_c));
% shift the signal to be centered mid frame
w = exp(-1i * 2 * pi * (n-1)/2 * [0:floor(n/2)-1 floor(-n/2):-1]'/ n); 


figure(3)
clf
subplot(411)
plot(real(Y(:,delay_index:i_max)))
title('azimuth real data at reflector')
legend

subplot(412)
plot(real(Ymf(:,delay_index:i_max)))
title('azimuth fast time match filtered real data at reflector')
legend

subplot(413)
plot(real(az))
title('expected real signal in Az')
legend
if mod(n, 2) == 0
	% Force conjugate symmetry. Otherwise this frequency component has no
	% corresponding negative frequency to cancel out its imaginary part.
	w(n/2+1) = real(w(n/2+1));
end 
Zmf = ifft( fft(z)  .* conj( fft(az) ) .* w);

subplot(414)
plot(abs(Zmf))
legend
% 
% figure(4)
% clf
% plot((R-depth)/lambda_c)
% title('Range change in units of wavelength')






