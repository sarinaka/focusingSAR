% Date last modified: 5/7/20
% Date created: 3/17/20

% script to run the radar simulator, should simulate the radar response
% from a point target, under different radar system conditions

clearvars;
close all

% user enabled/defined parameters
lambda = 1; % m
fc = 300e6; % Hz
tau_range = 50e-6; % sec
bw = 5e6; % Hz
slope_range = bw / tau_range; % Hz/sec
fs = 20e6; % Hz

vel_radar = 10; % m/s, radar platform velocity
d_target = 3000; % m, target depth
alpha = 20; % dB/km, attenuation rate
alpha = alpha / 1000; % linear/m

npts = 1000; % number of time steps in the simulation
theta_hbw = 30; % degrees, antenna 1/2 beamwidth angle
theta_hbw = theta_hbw * pi / 180;
sigma_noise = 0.1;

c = 3e8; % m/s

eps_ice = 3.17;
eps_bed = 80; % water
vel_ice = c / sqrt(eps_ice); % speed of wave in ice

% calculate reflection coefficient
r12 = (sqrt(eps_ice) - sqrt(eps_bed)) ./ (sqrt(eps_ice) + sqrt(eps_bed));
gamma = r12 .* conj(r12);

% min and max x distances for which the radar will illuminate the target
x_radar_min = -1 * tan(theta_hbw) * d_target;
x_radar_max = tan(theta_hbw) * d_target;

% time the target is illuminated
t_illum = linspace(x_radar_min/vel_radar, x_radar_max/vel_radar, npts);

x_radar = vel_radar * t_illum;
r_target = sqrt(x_radar.^2 + d_target^2); % range from radar to target, function of time

figure; plot(t_illum, r_target); xlabel('Time (s)'); ylabel('Distance (m)'); title('Range to Target')

%% make the reference chirp
N = 2^(ceil(log2(fs * tau_range))); % make reference chirp power of 2 in length for faster fft
N = N*2;
reference_chirp = makeChirp(slope_range, tau_range, fs, fc, N, 0);

t = linspace(0, tau_range, length(reference_chirp));

% take fft of reference chirp
reference_ft = fft(reference_chirp);
freq = (-fs/2:fs/N:(fs/2-fs/N)) + fc; % center frequency array at fc

figure
subplot(2,1,1)
hold on
plot(t, real(reference_chirp));
plot(t, imag(reference_chirp));
hold off
xlabel('Time (seconds)')
title('Transmitted Signal')
legend('Real Component', 'Imaginary Component')

subplot(2,1,2)
plot(freq/1e6, 20*log10(fftshift(abs(reference_ft))))
xlabel('Frequency (MHz)')
ylabel('Power (dB)')

%% make received chirps (delay, sum multiple copies, add noise, attenuate)

% calculate delay based on range vector
%delay = 2 * r_target / vel_ice; % two way travel time in ice
rx_signal = makeChirp(slope_range, tau_range, fs, fc, N, min(r_target));
rx_signal = rx_signal + makeChirp(slope_range, tau_range, fs, fc, N, 0); % add surface reflection

% add noise
rx_signal = rx_signal + sigma_noise * rand(size(rx_signal));

% vector that stores the attenuation in linear units to each range (kind of
% sketchy implemenation
attenuation = 10^(-25/20); % 25 dB/km attenuation, not sure whether it should be - or +

% apply attenuation
rx_signal = rx_signal * exp(-2 * attenuation / 1000 * min(r_target));

rx_signal_ft = fft(rx_signal);
pc_signal_ft = rx_signal_ft .* conj(reference_ft); % pulse compress signal
pc_signal = ifft(pc_signal_ft);

t_pc = linspace(-tau_range/2, tau_range/2, length(reference_chirp));

figure
subplot(2,1,1)
hold on
plot(t, real(reference_chirp));
plot(t, real(rx_signal));
hold off
xlabel('Time (seconds)')
title('Delayed Chirp Signal')
legend('Transmitted Signal', 'Received Signal')
% TODO: I think a lot of my axes labels are wrong in these plots
subplot(2,1,2)
plot(t_pc*1e6, 20*log10(abs(ifftshift(pc_signal))))
xlabel('Time (\mus)')
ylabel('Power (dB)')
title('Pulse Compressed Signal')

%% now create a record of the chirp for all range points

 % vector that stores the attenuation in linear units to each range (kind of
 % sketchy implemenation
 attenuation = 10^(-25/20); % 25 dB/km attenuation, not sure whether it should be - or +

% start with the dumb way of doing it--loop through all range bins and
% stack results in matrix (there's smarter way to do it with vectorization I'm
% sure)
rx_full_pc = zeros(N, length(r_target));
for ii = 1:length(r_target)
    rx_sig = makeChirp(slope_range, tau_range, fs, fc, N, r_target(ii));
    rx_sig = rx_sig + makeChirp(slope_range, tau_range, fs, fc, N, 0); % add surface reflection
    
    % add noise
    rx_sig = rx_sig + sigma_noise * rand(size(rx_sig));
    
    % apply attenuation
    rx_sig = rx_sig * exp(-2 * attenuation / 1000 * r_target(ii));
    
    % pulse compress
    rx_sig_ft = fft(rx_sig);
    pc_sig_ft = rx_sig_ft .* conj(reference_ft); % pulse compress signal
    pc_sig = ifft(pc_sig_ft);
    
    % add to data array
    rx_full_pc(:,ii) = pc_sig; 
end

% TODO: y axis is offby a factor of 2, not sure whether it's a label issue
% or a real missing 2 somewhere
% display fake range compressed radargram
figure
imagesc(x_radar, t_pc*vel_ice/2, 20*log10(abs(ifftshift(rx_full_pc,1))))
xlabel('Along Track Distance Relative to Scatterer (m)')
ylabel('Depth (m)')
c = colorbar;
c.Label.String = 'Power (dB)';


%% do range migration


%% do azimuth compression

