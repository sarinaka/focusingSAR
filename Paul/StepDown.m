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
f_c = 10e6;           % Center frequency of chirp [Hz] 
BW = 20e6;           % BW of chirp [Hz]
f_0 = f_c - BW/2;    % Initial frequency [Hz]           
t_c = 6e-5;           % Chirp Length [s]
% Survey Parameters
c = 3e8/1.79;
lambda_c = c/f_c;     % Wavelength (center) in ice [m]
lambda_0 = c/f_0;     % Wavelength (initial) in ice [m]
depth = 4e3;          % Scatter Depth [m]
offset = 1e3;
n = 101;
dx = offset/n;
xx = -dx*(n-1)/2:dx:dx*(n-1)/2;

% Plotting params
delay_time = 2*depth/c;

%% Make pulse
t_sub = 0:1/f_s:t_c; %pulse only for duration of pulse
X = zeros(size(t));  %pulse signal is 0 otherwise
X(1:length(t_sub)) = exp(1i*(pi.*(BW/t_c).*t_sub.^2+2.*pi.*f_0.*t_sub)); %LFMCM pulse

disp("range resolution is " + c/(2*BW) + " m");
disp("deramped freqency is " + 2*BW/(t_c*c) + "Hz/meter => " +2*depth*BW/(t_c*c)/1e6 + " MHz" )

%Prep fields for output, range (L) is row, azimuth (n) is column)
%BACKWARDS OF De Wit 2016. I also plot everything as an ' which is 
Y = zeros(L,n);
Ysd = zeros(L,n);

%Loop over surface points
for i = 1:n
    %Find range, shift chirp (in Time and Fx)
    r1 = sqrt((depth).^2 + ((i-(n+1)/2)*dx).^2);
%     r2 = sqrt((depth+4e1).^2 + ((i-(n+1)/(2.5))*dx).^2);
    y_tmp = chirpOut(X,t,r1,0,f_c,f_s);% + chirpOut(X,t,r2,0,f_c,f_s);
    
    %raw data
    Y(:,i) = y_tmp;
    Ysd(:,i) = real(y_tmp).*X;
    clear r y_tmp
end

I = processBlock(Y.',X,f_s,f_c,dx,t);   %<--- watch out to not do complex cong transpose





MF = ifft( fft(Y(:,(n+1)/2))  .* conj( fft(X) ) ); %should scale ref to 1


%% FFT in range
% x = fftshift(X);
% y = fftshift(Y(:,(n+1)/2));
% ysd = fftshift(Ysd);
% f = f_s*(-(L-1)/2:(L-1)/2)/L/1e6;

x = fft(X);
y = fft(Y(:,(n+1)/2));
ysd = fft(Ysd);
%clip for posf requency values only
ysdclp = ysd(1:(L+1)/2,:);
f = f_s*([0:(L-1)/2,-(L-1)/2:-1])/L/1e6;
f_plus = f_s*(0:(L-1)/2)/L;

% Range migrate up for items on center frame
[rg,az] = ndgrid(f_plus/(2*BW/(t_c*c)),xx);
Rshift = griddedInterpolant(rg,az,ysdclp,'linear','none');
% Interpolate signal back up
ysdrm = Rshift(sqrt((f_plus./(2*BW/(t_c*c)))'.^2 + az.^2),az); %exact

% Az focus



%% 
figure(2)
clf
subplot(411)
    plot(t,real(X))
    legend('Tx')
    xlabel('Time')
subplot(412)
    plot(t,real(Y(:,(n+1)/2)))
    legend('Rx')
    xlabel('Time')
subplot(413)
    plot(f,abs(x),'-*')
    hold on
    plot(f,abs(y))
    plot(f,abs(ysd(:,(n+1)/2)))
    legend('Tx','Rx','SD')
    xlabel('Frequency [MHz]')
subplot(414)
    plot(t,real(X))
    hold on
    plot(t,abs(MF)/max(abs(MF)))
    plot(t,real(Y(:,(n+1)/2)))
    plot(t,real(Ysd(:,(n+1)/2)))
    legend('Tx','MF','Rx','SD')
    xlabel('Time')

figure(3)
    subplot(411)
        prettyPlot(xx,t*c/2,real(Y))
%         ylim([delay_time-.1e-5 delay_time+.1e-5])
        ylabel('range 1 way [m]')
        xlabel('along track')
        title('Raw data')
        colorbar
    subplot(412)
        prettyPlot(xx,t*c/2,abs(I'))
%         ylim([delay_time-.1e-5 delay_time+.1e-5])
        colorbar
        ylabel('range 1 way [m]')
        xlabel('along track')
        title('Full Focused Data')
    subplot(413)
        prettyPlot(xx,f_plus/(2*BW/(t_c*c)),abs(ysdclp))
%         ylim([delay_time-.1e-5 delay_time+.1e-5])
        colorbar
        ylabel('range [m]')
        xlabel('along track')
        title('Stepped Down')
    subplot(414)
        prettyPlot(xx,f_plus/(2*BW/(t_c*c)),abs(ysdrm))
        ylim([depth-50,depth+50])
        colorbar
        ylabel('range [m]')
        xlabel('along track')
        title('Stepped Down migrated')
