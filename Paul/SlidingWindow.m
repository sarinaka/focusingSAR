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
BW = 20e6;           % BW of chirp [Hz]
f_0 = f_c - BW/2;    % Intial frequency [Hz]           
t_c = 10e-5;           % Chirp Length [s]
% Survey Parameters
c = 3e8/1.31;
lambda_c = c/f_c;     % Wavelength (center) in ice [m]
lambda_0 = c/f_0;     % Wavelength (initial) in ice [m]
n = 251; 			  % surface sample points
dx = 5;  	          % Distance between sample points [m]
xx = (((1:n)*dx)-(n+1)/2*dx)';
depth = 4e3;          % Scatter Depth [m]
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
    r2 = sqrt((depth+4e1).^2 + ((i-(n+1)/(2.5))*dx).^2);
    y_tmp = chirpOut(X,t,r1,0,f_c,f_s) + chirpOut(X,t,r2,0,f_c,f_s);
    
    %raw data
    Y(i,:) = y_tmp; 
    Ysd(i,:) = y_tmp.*X;
    clear r y_tmp
end
clear i

tic

I = zeros(size(Y));
Isd = zeros(size(Y));
for i = 1:1:((n-ap)+1)
    tmp = processBlock(Y(i:(ap+i-1),:),X,f_s,f_c,dx,t)./window; %Processing Az with initial freq (not center) works way better
    I(((ap-window)/2+i):(ap+i-1-(ap-window)/2),:) = ...
        I(((ap-window)/2+i):(ap+i-1-(ap-window)/2),:) + tmp((ap-window)/2+1:end-(ap-window)/2,:);
    tmp = processBlock(Ysd(i:(ap+i-1),:),X,f_s,f_c,dx,t)./window; %Processing Az with initial freq (not center) works way better
    Isd(((ap-window)/2+i):(ap+i-1-(ap-window)/2),:) = ...
        Isd(((ap-window)/2+i):(ap+i-1-(ap-window)/2),:) + tmp((ap-window)/2+1:end-(ap-window)/2,:);
end
toc
%% Plot the returns, abs() of complex values to display
figure(1)
clf
subplot(221)
	prettyPlot(real(Y'))
	ylabel('range')
	xlabel('along track')
	title('Raw data')
	colorbar
subplot(223)
    prettyPlot(abs(I(:,i_min:i_max)'))
    colorbar
    ylabel('range')
    xlabel('along track')
    title('Full Focused Data')
subplot(222)
	prettyPlot(real(Ysd'))
	ylabel('range')
	xlabel('along track')
	title('Raw data')
	colorbar
subplot(224)
    prettyPlot(abs(Isd(:,i_min:i_max)'))
    colorbar
    ylabel('range')
    xlabel('along track')
    title('Full Focused Data')
    
disp("SNR is " + round(10*log10(max(max(abs(I)))/mean(mean(abs(I)),'omitnan')),2) + "dB");
    
% figure(2)
% clf
% subplot(211)
%     plot(t,real(X)/max(abs(X)))
%     hold on
%     plot(t,real(Y((n+1)/2,:))/max(abs(Y((n+1)/2,:))))
%     legend('Chirp', 'Raw data')
% subplot(212)
%     plot(t,abs(I((n+1)/2,:))/max(abs(I((n+1)/2,:))))
%     hold on
%     plot(t,real(X.*Y((n+1)/2,:)));
%     legend('Focused','Stepped Down')
%     sgtitle('Center point fast time down')

