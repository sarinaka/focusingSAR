function vdat = LongBurstRMB3(Filename, Burst, SamplesPerChirp, MaxNChirp)

% vdat = LongBurstRMB3(Filename, Burst, SamplesPerChirp, MaxNChirp)
%
% Load data from RMB-A and average by stacking to reduce memory use
% Used for averaging very long bursts

WperChirpHdr = 0;
MaxHeaderLen = 500;
burstpointer = 0;
vdat.Code = 0;
NChirp = 0;
vdat.v = zeros(1,SamplesPerChirp);
fid = fopen(Filename,'r');
if fid >= 0
    fseek(fid,0,'eof');
    filelength = ftell(fid);
    BurstCount = 1;
    while BurstCount <= Burst && burstpointer <= filelength - MaxHeaderLen
        fseek(fid,burstpointer,'bof');
        A = fread(fid,MaxHeaderLen,'*char');
        A = A';
        searchind = strfind(A, 'Samples:');
        if ~isempty(searchind)
            try
                searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
                vdat.Nsamples = sscanf(A(searchind(1)+8:searchCR(1)+searchind(1)),'%d');
                WperChirpCycle = vdat.Nsamples + WperChirpHdr;
                searchind = strfind(A, 'Chirps in burst:');
                searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
                vdat.ChirpsInBurst = sscanf(A(searchind(1)+16:searchCR(1)+searchind(1)),'%d');
                
                searchind = strfind(A, '*** End Header ***');
                
                burstpointer = burstpointer + searchind(1) + 20;
            catch
                vdat.Code = -2;
                vdat.Burst = BurstCount;
                return
            end
        end
        WordsPerBurst = vdat.ChirpsInBurst * WperChirpCycle;
        if BurstCount < Burst && burstpointer <= filelength - MaxHeaderLen
            burstpointer = burstpointer + vdat.ChirpsInBurst * WperChirpCycle*2;
        end
        BurstCount = BurstCount + 1;
    end
    
    % Extract remaining information from header
    searchind = strfind(A, 'Time stamp:');
    searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
    try
        td = sscanf(A(searchind(1)+11:searchCR(1)+searchind(1)),...
            '%d-%d-%d %d:%d:%d');
        vdat.TimeStamp = datenum(td(1),td(2),td(3),td(4),td(5),td(6));
    catch err
        vdat.Code = 1;
    end
    
    searchind = strfind(A, 'Temperature 1:');
    try
        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        vdat.Temperature_1 = sscanf(A(searchind(1)+14:searchCR(1)+searchind(1)),'%f');
    catch err
        vdat.Code = 1;
    end
    
    searchind = strfind(A, 'Temperature 2:');
    try
        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        vdat.Temperature_2 = sscanf(A(searchind(1)+14:searchCR(1)+searchind(1)),'%f');
    catch err
        vdat.Code = 1;
    end
    
    searchind = strfind(A, 'Battery voltage:');
    try
        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        vdat.BatteryVoltage = sscanf(A(searchind(1)+16:searchCR(1)+searchind(1)),'%f');
    catch err
        vdat.Code = 1;
    end
    
    searchind = strfind(A, 'Attenuator 1:');
    try
        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        vdat.Attenuator_1 = sscanf(A(searchind(1)+13:searchCR(1)+searchind(1)),'%f',4);
    catch err
        vdat.Code = 1;
    end
    
    searchind = strfind(A, 'Attenuator 2:');
    try
        searchCR = strfind(A(searchind(1):end),[char(13),char(10)]);
        vdat.Attenuator_2 = sscanf(A(searchind(1)+13:searchCR(1)+searchind(1)),'%f',4);
    catch err
        vdat.Code = 1;
    end
    
    fseek(fid,burstpointer-1,'bof');
    if BurstCount == Burst+1
        while vdat.Code == 0 && NChirp < vdat.ChirpsInBurst && NChirp < MaxNChirp
            [tmp count] = fread(fid,SamplesPerChirp,'*int16','ieee-le');
            if count < SamplesPerChirp
                vdat.Code = 2;
            else
                vdat.v = vdat.v + double(tmp');
                NChirp = NChirp + 1;
                fseek(fid,burstpointer-1+vdat.Nsamples * NChirp * 2,'bof');
            end
        end
        %vdat.v = single(vdat.v/NChirp);
        vdat.v = vdat.v * 2.5 / 2^16 + 1.25;
        vdat.Startind = 1;
        vdat.Endind = vdat.Startind + SamplesPerChirp - 1;
        vdat.Burst = Burst;
        vdat.NChirp = NChirp; % save number of chirps so we can average multiple bursts
    else
% Too few bursts in file
        vdat.Burst = BurstCount - 1;
        vdat.Code = -4;
    end
    fclose(fid);
else
% Unknown file
    vdat.Code = -1;
end
