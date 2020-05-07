% This is a script to make and process n chirps from stationary platform, while moving 
% between chirps to make a scene. Produces both raw data and match filtered
% data products. The scene is focused with a sliding window that is the
% synthetic aperture, which is less than the total survey length.
clear

%% Inputs
c = 3e8/1.79;         % Speed of light in ice [m/s]
% Rx parameters
f_s = 60e6;           % Sampling frequency [Hz]
T   = 10e-5;          % Record Time [s]
t = 0:1/f_s:T;        % Time vector [s]
L = length(t);        % Recording vector length [ ]
% Tx parameters
f_c = 10e6;           % Center frequency of chirp [Hz] 
BW = 20e6;            % BW of chirp [Hz]
f_0 = f_c - BW/2;     % Initial frequency [Hz]           
t_c = 5e-5;           % Chirp Length [s]
% Survey parameters
depth = 4e3;          % Scatter Depth [m]
theta = pi/12;        % beam half width [rad]
surveyVelocity = 10;  % Velocity of platform [m/s](still assuming point and shoot)
PRF = 2;              % pulse repetition frequency [Hz]
f_zone = 2*tan(theta)*depth; %distance scatter is visible
% Synthetic data parameters
dx = surveyVelocity/PRF; % Dist between sample points [m]
n = floor(f_zone/dx/2)*2+1; 	 % Surface sample points (must be odd) 	         
xx = (((1:n)*dx)-(n+1)/2*dx)'; % sample points
% Focusing Parameters
SynAp = 400;                 % Aperture width [m]
ap = floor(SynAp/dx/2)*2+1;  % Aperture width in Az bins (must be odd)
window = ap;                 % Window width in Az bins
step = 1;                    % stepping size for plotting processed looks 
% Plotting parameters
delay_time = 2*depth/c;

disp(round(f_zone/1e3,2) + "km survey with " + n + " samples");
%% Make pulse
t_sub = 0:1/f_s:t_c; %pulse only for duration of pulse
X = zeros(size(t));  %pulse signal is 0 otherwise
X(1:length(t_sub)) = exp(1i*(pi.*(BW/t_c).*t_sub.^2+2.*pi.*f_0.*t_sub)); %LFMCM pulse

%% Sweep along surface, let full pulse out and wait for return, then move.
Y = zeros(n,L); %initialize raw data field
Ysd = Y;        %initialize step down data field
R = zeros(n,1);
for i = 1:n
    %Find range, shift chirp (in Time and Fx)
    r1 = sqrt((depth).^2 + ((i-(n+1)/2)*dx).^2);
%     r2 = sqrt((depth+4e1).^2 + ((i-(n+1)/(2.5))*dx).^2);
    y_tmp = chirpOut(X,t,r1,0,f_c,f_s);% + chirpOut(X,t,r2,0,f_c,f_s);
    
    %raw data
    Y(i,:) = y_tmp;
    clear r y_tmp
end
clear i

tic

Imf = zeros(size(Y));
Ish = zeros(size(Y));
I = zeros(size(Y));
for i = 1:step:((n-ap)+1)
    [tmf,tmfsh,tmfshaz] = processBlock(Y(i:(ap+i-1),:),X,f_s,f_c,dx,t); 
    Imf(((ap-window)/2+i):(ap+i-1-(ap-window)/2),:) = ...
        Imf(((ap-window)/2+i):(ap+i-1-(ap-window)/2),:) + tmf((ap-window)/2+1:end-(ap-window)/2,:);
    Ish(((ap-window)/2+i):(ap+i-1-(ap-window)/2),:) = ...
        Ish(((ap-window)/2+i):(ap+i-1-(ap-window)/2),:) + tmfsh((ap-window)/2+1:end-(ap-window)/2,:);
    I(((ap-window)/2+i):(ap+i-1-(ap-window)/2),:) = ...
        I(((ap-window)/2+i):(ap+i-1-(ap-window)/2),:) + tmfshaz((ap-window)/2+1:end-(ap-window)/2,:);
end
toc
%% Plot the returns, abs() of complex values to display
figure(1)
set(gcf,'Position',[100 100 600 1200])
clf
subplot(411)
	prettyPlot(xx,t,real(Y'))
	ylabel('range')
	xlabel('along track')
	title('Raw data')
	colorbar
subplot(412)
    prettyPlot(xx,t,abs(Imf'))
    ylim([delay_time-.1e-5 delay_time+.1e-5])
    colorbar
    ylabel('range')
    xlabel('along track')
    title('Match Filter Data')
subplot(413)
    prettyPlot(xx,t,abs(Ish'))
    ylim([delay_time-.1e-5 delay_time+.1e-5])
    colorbar
    ylabel('range')
    xlabel('along track')
    title('Range Shifted Data')
subplot(414)
    prettyPlot(xx,t,abs(I'))
    ylim([delay_time-.1e-5 delay_time+.1e-5])
    colorbar
    ylabel('range')
    xlabel('along track')
    title('Full Focused Data')
sgtitle("F_s = " + f_s/1e6 + " MHz, F_c = " + f_c/1e6 + " MHz, BW = " + BW/1e6 + " MHz")    
disp("SNR is " + round(10*log10(max(max(abs(I)))/mean(mean(abs(I)),'omitnan')),2) + "dB");
    
figure(2)
clf
subplot(211)
    plot(t,real(X))
    hold on
    plot(t,real(Y((n+1)/2,:)))
    legend('Chirp', 'Raw data')
subplot(212)
    plot(t,abs(I((n+1)/2,:))/max(abs(I((n+1)/2,:))))
    hold on
    plot(t,real(X.*Y((n+1)/2,:)));
    legend('Focused','Stepped Down')
    sgtitle('Center point fast time down')
