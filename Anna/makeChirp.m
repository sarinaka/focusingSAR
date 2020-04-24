function [chirp_final] = makeChirp(slope, tau, fs, fc, N, delay)
% makeChirp creates an FMCW chirp with given parameters
% INPUTS: 
    % slope: chirp slope, Hz/sec
    % tau: pulse length, sec
    % fs: sample frequency, Hz
    % fc: chirp center frequency, Hz
    % N: total number of samples in chirp (length including any zeros)
    % delay: range to delay the chirp by 
% OUTPUTS:
    % chirpFinal: complex chirp
% NOTES:
    % chirp bandwidth = slope * tau
    
c = 3e8; % m/s
eps_ice = 3.17; 
vel_ice = c / sqrt(eps_ice); % speed of wave in ice

npts = tau * fs; % number of points in the actual transmitted chirp
t = linspace(0, tau, npts); % time vector
%t = t - delay; 

f1 = fs - (slope * tau / 2); % get start frequency

chirp_phase = pi*slope*t.^2 + 2*pi*fc*t; 
chirp_data = exp(1i * chirp_phase);

%chirp_delayed = fft(chirp_data) * exp(-1j * 2 * delay / vel_ice);
%chirp_final = ifft(chirp_delayed);
delay = floor(2 * delay / vel_ice * fs);
chirp_final = [zeros(1, delay) chirp_data]; % delay chirp if desired
%chirp_final = chirp_data;
chirp_final(end+1:N) = 0; % zero pad to end of array

end

