function envFilter = computeTimbreFilter(targetTimbreEnv, inputSpec, inputPitch, params)
% function envFilter = computeTimbreFilter(targetTimbreEnv, inputSpec,
% inputPitch, params)
%
% Given the target timbre envelope and an input spectrum, it computes the
% necessary filter to preserve the timbre. It can be used for morphing two sounds.
%
    numMaxHarmonics = params.numMaxHarmonics ; %500;
    srate = params.fs;
    thresholdFilterdB = 20; % limit gain in dB
    thresF = 10^(thresholdFilterdB/20);
    crossFreq = 7.5e3;
    specSize = size(targetTimbreEnv,2);
    crossFreqBin = floor(crossFreq * 2 * size(targetTimbreEnv,2)/srate);
    rampBins = min(crossFreqBin, specSize- crossFreqBin);
    crossBand = [zeros(1,crossFreqBin), 1/rampBins:1/rampBins:1, ones(1,max(0,specSize - (crossFreqBin+rampBins))) ]; 
    crossBand = sqrt(crossBand);

    % compute envelope
    [synmagfftEnv_PS, ~,  ~, ~] = harmonicEnv(inputSpec , inputPitch, srate, srate/2, numMaxHarmonics);
    
    % smooth envelope in frequency above half the cross frequency   FIR
    % with h = (1/3)[1,1,1]
    targetTimbreEnv(:,crossFreqBin:specSize) = (1/4) * (targetTimbreEnv(:,crossFreqBin-3:specSize-3) + targetTimbreEnv(:,crossFreqBin-2:specSize-2) + ...
        targetTimbreEnv(:,crossFreqBin-1:specSize-1) + targetTimbreEnv(:,crossFreqBin:specSize));
    synmagfftEnv_PS(:,crossFreqBin:specSize) =  (1/4) * (synmagfftEnv_PS(:,crossFreqBin-3:specSize-3) +synmagfftEnv_PS(:,crossFreqBin-2:specSize-2) + ...
        synmagfftEnv_PS(:,crossFreqBin-1:specSize-1) + synmagfftEnv_PS(:,crossFreqBin:specSize));
    
    % compute filter
    envFilter = ones(size(targetTimbreEnv));
    validIdx = find(sum(synmagfftEnv_PS,2)>0);
    envFilter(validIdx,:) = (targetTimbreEnv(validIdx,:) ./ synmagfftEnv_PS(validIdx,:)); % get filter
    % high-shelf by-pass filter in order to preserve timbre only in low-mid
    % frequencies to avoid artifacts.
    
    %envFilter = envFilter .* (ones(size(envFilter,1),1) *  (1 - crossBand)) + ones(size(envFilter)).* (ones(size(envFilter,1),1) * crossBand); % commented test edwin
    envFilter = min(thresF,envFilter);
    % normalize energy
    eIn = sum(inputSpec.^2,2);    
    eFilt = sum((inputSpec.*envFilter).^2,2);
    scaleFact = min(thresF,sqrt(eIn./eFilt))*ones(1,size(envFilter,2));
    envFilter = envFilter .* scaleFact;
    % cross facde not to increase HF beyond the corssfade
    envFilter = envFilter .* (ones(size(envFilter,1),1) *  (1 - crossBand)) + ones(size(envFilter)).* (ones(size(envFilter,1),1) * crossBand); % commented test edwin
    
    
    % smooth over time to avoid rapid changes in the timbre filter
    h = [0.1 0.25 1 0.25 0.1];       
    %original function
    %mask = filter2(h/sum(h),envFilter')'; % smooth
    % alternative reduce memory requirement
    specSize = size(envFilter,2);
    specStep  = floor(specSize/128);    
    for j=1:specStep:specSize-specStep,        
        envFilter(:,j:j+specStep) = filter2(h/sum(h),envFilter(:,j:j+specStep)')';
    end
    
end
