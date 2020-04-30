% This is a script to make n chirps from stationary platform, while moving 
% between chirps to make a scene. Produces both raw data and match filtered
% data products
clear
%% Inputs
% Rx parameters
f_s = 60e6;          % Sampling frequency [Hz]
T   = 15e-5;           % Record Time [s];
t = 0:1/f_s:T;        % Time vector [s]
L = length(t);        % Recording vector length [ ]
% Tx parameters
f_c = 40e6;           % Center frequency of chirp [Hz] 
BW = 10e6;           % BW of chirp [Hz]
f_0 = f_c - BW/2;    % Intial frequency [Hz]           
t_c = 10e-5;           % Chirp Length [s]
% Survey Parameters
c = 3e8/1.31;
lambda_c = c/f_c;     % Wavelength (center) in ice [m]
lambda_0 = c/f_0;     % Wavelength (initial) in ice [m]
n = 5; 			  % surface sample points
dx = 5;  	          % Distance between sample points [m]
xx = (((1:n)*dx)-(n+1)/2*dx)';
depth = 1e3;          % Scatter Depth [m]
% Focusing Parameters
ap = 51;              % Apeture width in Az bins
window = 5;           % Window width in Az bins

% Plotting params
delay_index = round(depth*2/c*f_s);
i_range = 100;
i_min = delay_index-i_range;
i_max = delay_index+i_range;

%% Make pulse
t_sub = 0:1/f_s:t_c; %pulse only for duration of pulse
X = zeros(size(t));  %pulse signal is 0 otherwise
X(1:length(t_sub)) = exp(1i*(pi.*(BW/t_c).*t_sub.^2+2.*pi.*f_0.*t_sub)); %LFMCM pulse

%% Sweep along surface, let full pulse out and wait for return, then move.
Y = zeros(n,L); %initialize raw data field
Ysd = Y;        %initialize stepdown data field
R = zeros(n,1);
for i = 1:n
    %Find range, shift chirp (in Time and Fx)
    r1 = sqrt((depth).^2 + ((i-(n+1)/2)*dx).^2);
    y_tmp = chirpOut(X,t,r1,0,f_c,f_s);
    
    %raw data
    Y(i,:) = y_tmp; 
    Ysd(i,:) = y_tmp.*X;
    clear r y_tmp
end
y = fft(Y(3,:));
ysd = fft(Ysd(3,:));
f = f_s*(0:(L-1))/L;

%% 
figure(1)
clf
subplot(211)
    plot(f,abs(y(:)))
    hold on
    plot(f,abs(ysd(:)))
    legend('Tx','SD')
subplot(212)
    plot(t,real(X))
    hold on
    plot(t,real(Ysd(3,:)))
    legend('Tx','SD')

