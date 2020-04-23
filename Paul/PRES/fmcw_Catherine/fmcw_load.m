function vdat = fmcw_load(filename,permittivity,burst)

% vdat = fmcw_load(filename,burst)
%
% Load FMCW radar burst and metadata
%
% input: filename and burst number
% output: vdat - structure of metadata and data

% Craig Stewart
% 2013 April 24
% 2013 September 30 - corrected error in vif scaling
% 2014/5/20 time stamp moved here from fmcw_derive_parameters (so that this
% is not overwritted ater)
% 2014/5/21 changed howl radar chirp is defined (now using chirp gradient as
% fundamental parameter)
% 2014/5/22 fixed bug in chirptime
% 2014/8/21 moved make odd length to external (called from fmcw_range)
% 2014/10/22 KWN - edited to allow for new headers in RMB2 files

%permittivity = 2.9;
% Check args
if nargin == 0
    [filename, path] = uigetfile(['*.dat;*.DAT;*.000;*.mat'],'Choose radar file','multiselect','off');
    if isa(filename,'double') % no files chosen
        return
    end
    filename = [path,filename];
end
if nargin < 3
    burst = 1;
end

% Settings
SamplingFrequency = 40000; %
doShowFileType = 0;

%% Load data and reshape array
[path,~,ext] = fileparts(filename);
if strcmp(ext,'.mat')
    load(filename); % note when you load a mat file you get whatever burst was stored in this - not the one you selected
    FileFormat = 'mat';
else
    FileFormat = fmcw_file_format(filename);
    if FileFormat == 5
        vdat = LoadBurstRMB5(filename, burst, SamplingFrequency);% Data from after Oct 2014 (RMB2b + VAB Iss C, SW Issue >= 101)
    elseif FileFormat == 4
        vdat = LoadBurstRMB4(filename, burst, SamplingFrequency); % Data from after Oct 2013  (RMB1b)
        %disp('file format4')
    elseif FileFormat == 3
        vdat = LoadBurstRMB3(filename, burst, SamplingFrequency); % Data from Jan 2013 (RMB1a)
    elseif FileFormat == 2
        %vdat = LoadOldBurstCompat(filename, burst, SamplingFrequency); % Data from Prototype FMCW radar (nov 2012) (RMBA?)
        vdat = LoadBurstMcM(filename,burst);
    end
end
vdat.FileFormat = FileFormat;
if doShowFileType
    disp(['Detected file format: ' int2str(FileFormat)])
end

% Check file was found
switch(vdat.Code)
    case -1
        fprintf('Unable to open file: %s\n',filename);
        return
    case -2
        fprintf('Corrupt header in burst - fatal %d\n',vdat.Burst);
        return
    case -4
        %        disp(['Burst ' int2str(burst) ' not found in file ' filename]);
        return
end

% Extract just good chirp data from voltage record and rearrange into
% matrix with one chirp per row
% note: you can't just use reshape as we are also cropping the 20K samples
% of sync tone etc which occur after each 40K of chirp.
AttSet = vdat.Attenuator_1 + 1i*vdat.Attenuator_2; % unique code for attenuator setting


%% Add metadata to structure
% Sampling parameters
vdat.filename = filename;
if ~ischar(FileFormat)
    vdat.SamplesPerChirp = vdat.Nsamples;
    vdat.fs = 4e4; % sampling frequency
    vdat.f0 = 2e8; % start frequency
    %vdat.fc = 3e8; % start frequency
    vdat.K = 2*pi*2e8; % chirp gradient in rad/s/s (200MHz/s)
    %vdat.f0 = vdat.f0 + (vdat.K/(4*pi))/vdat.fs; % start frequency
    vdat.processing = {};
    
    if FileFormat == 5
        H = fmcw_ParametersRMB2(vdat.filename);
    elseif FileFormat == 4
        H = fmcw_ParametersRMB1b(vdat.filename);
    end
    if FileFormat == 5 || FileFormat == 4
        vdat.K = H.K;
        vdat.f0 = H.startFreq;
        vdat.fs = H.fs;
        vdat.f1 = H.startFreq + H.chirpLength * H.K/2/pi;
        vdat.SamplesPerChirp = round(H.chirpLength * H.fs);
        vdat.T = H.chirpLength;
        vdat.B = H.chirpLength * H.K/2/pi;
        vdat.fc = H.startFreq + vdat.B/2;
        vdat.dt = 1/H.fs;
        vdat.er = permittivity;
%         vdat.er = 3.1;
        vdat.ci = 3e8/sqrt(vdat.er);
        vdat.lambdac = vdat.ci/vdat.fc;
        vdat.Nsamples = H.nchirpSamples;
    else
        vdat.er = permittivity;
%         vdat.er = 3.1;
        vdat = fmcw_derive_parameters(vdat);
    end
else
    vdat.er = permittivity;
%     vdat.er = 3.1;
    vdat.dt = 1/vdat.fs;
    vdat.ci = 3e8/sqrt(vdat.er);
    vdat.lambdac = vdat.ci/vdat.fc;
end

vdat.er
vdat.Endind = vdat.Startind + vdat.SamplesPerChirp - 1;

% Load each chirp into a row
vdat.vif = zeros(vdat.ChirpsInBurst,vdat.SamplesPerChirp); % preallocate array
chirpInterval = 1.6384/(24*3600); % days
for chirp = 1:vdat.ChirpsInBurst
    vdat.vif(chirp,:) = vdat.v(vdat.Startind(chirp):vdat.Endind(chirp));
    vdat.chirpNum(chirp,1) = chirp; % chirp number in burst
    vdat.chirpAtt(chirp,1) = AttSet(1+mod(chirp-1,numel(AttSet))); % attenuator setting for chirp
    vdat.chirpTime(chirp,1) = vdat.TimeStamp + chirpInterval*(chirp-1); % time of chirp
end
vdat.er
% Create time and frequency stamp for samples
vdat.t = vdat.dt*(0:size(vdat.vif,2)-1); % sampling times (rel to first)
vdat.f = vdat.f0 + vdat.t.*vdat.K/(2*pi);
vdat.er
% Calibrate
%ca13 = [1 6]; % 2013
%ca14 = [1 2]; % 2014
%ca = [1 4];
%vdat = fmcw_cal(vdat,ca13);
