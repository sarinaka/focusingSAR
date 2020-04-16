% Make n chirps while moving btwn chirps, try making a radargram

clear
% Rx
F_s = 3e4;            % Sampling frequency [Hz]
T   = .1;              %Record Time [s];
t = 0:1/F_s:T;       % Time vector [s]
L = length(t);      % Recording length [ ]
% Tx
f_c = 1e2;            % Intial frequency of chirp [Hz] 
swp = 1.5e5;           % Sweep frequency of chirp [Hz/s]
t_c = 1e-2;            % Chirp Length [s]

n = 100; %surface sample points
Y = zeros(n,L);
Ymf = Y;

t_sub = 0:1/F_s:t_c;
X = zeros(size(t));
X(1:length(t_sub)) = exp(1i*(pi.*swp.*t_sub.^2+2.*pi.*f_c.*t_sub));

%sweep along surface, let full pulse out and wait for return, then move.
for i = 1:n
    %Find range, shift chirp
    r = sqrt((5e6).^2 + ((i-50)*2e5).^2);
    y_tmp = chirpOut(X,t,r,0,f_c);
    
    %raw data
    Y(i,:) = y_tmp; 
    
    % pulse compress in range, this could be done in a separate processing
    % loop if we wanted, indepentent to data collection
    Ymf(i,:) = ifft( fft(y_tmp)  .* conj( fft(X) ) );
end

%% Plot the returns, abs() of complex values to display
% figure(1) %Faster plotting option, but can't see waves
% clf
% subplot(211)
% imagesc(abs(Y'))
% ylabel('range')
% xlabel('along track')
% title('Raw data')
% colorbar
% subplot(212)
% imagesc(abs(Ymf'))
% ylabel('range')
% xlabel('along track')
% title('Match Filtered data')
% colorbar

figure(1) %very slow plotting, but cool to see waveforms!
clf
subplot(211)
wiggle(real(Y'))
ylabel('range')
xlabel('along track')
title('Raw data')
subplot(212)
wiggle(abs(Ymf'))
ylabel('range')
xlabel('along track')
title('Match Filtered data')