function [outResMag outResPhase] = computeResidualFromHarmMask(magfft, transmagfft, transphase, pitch, harmonicMask, nHarmonics, params)
% function [outResMag outResPhase] = computeResidualFromHarmMask(magfft, transmagfft, transphase, pitch, harmonicMask, nHarmonics, params)
%
%
% Computes the residual from an harmonic mask


specsize = size(magfft,2);
freq2bin = 2*specsize/params.fs;
nFr = size(magfft,1);
maxFreqBin = specsize;

outResMag = zeros(size(magfft));
outResPhase = zeros(size(magfft));
for  i=1:nFr
    
    if pitch(i)> 0 % 60Hz
        
        % shifted harmonics to find valleys between partials
        harmIdx = 1 + round([pitch(i)/2: pitch(i):params.fs/2]*freq2bin); % 0Hz correposnd to 1 for matlalb indexing
        harmIdx = harmIdx(1:min(length(harmIdx),nHarmonics));
        lastEnvBin = harmIdx(end);
        magDB = 20*log10(magfft(i,:));
        itype = 'linear';
        f = interp1([1 harmIdx maxFreqBin],[magDB(harmIdx(1)), magDB(harmIdx), min(magDB(harmIdx))],1:lastEnvBin,itype);
        flin = 10.^(f/20);
        outResMag(i,1:lastEnvBin) = [flin.*(flin>0)];
        
        outResMag(i,lastEnvBin+1:maxFreqBin) = transmagfft(i,lastEnvBin+1:maxFreqBin);
        
        outResPhase(i,1:lastEnvBin) = transphase(i,1:lastEnvBin) .* (1- harmonicMask(i,1:lastEnvBin)) +  harmonicMask(i,1:lastEnvBin) .* (2*pi*(-0.5 + rand(1,lastEnvBin)));
        outResPhase(i,lastEnvBin+1:maxFreqBin) = transphase(i,lastEnvBin+1:maxFreqBin);
    else
        outResMag(i,:) = transmagfft(i,:);
        outResPhase(i,:) = transphase(i,:);
    end
end
