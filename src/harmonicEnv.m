% compute the harmonic envelope for all frames given a magnitude
% spectrogram and frame pitch values
function [magfftEnv excitSlope harmAmps, harmIdxs]= harmonicEnv(magfft, pitch, Fs, forceExSlopeFreq, numMaxHarmonics)
% - magfft: matrix containing the magnitude spectrum of all frames (rows: num frames, cols: specsize)
% - pitch: column vector containing the pitch value for all frames
% - Fs: sampling rate 
% - forceExcitSlopeFr: freq above which the residual slope is forced to an exponential decay (in Hz),
%                     if set to -1, the original spectral envelope is
%                     returned. (default 5kHz). 
% - numMaxHarmonics: number of maximuim aronic partials detected.
% Outputs:
% - magfftEnv:  matrix containing the magnitude envelope spectrum of all frames
% (rows: num frames, cols: specsize)
% - excSlope: vector estimated excitation slope in dB/octave per frame

%minPitchHz = 35; % Hz
specsize = size(magfft,2);
freq2bin = 2*specsize/Fs;
maxFreqBin = specsize; 
magfftEnv = zeros(size(magfft,1), maxFreqBin);
excitSlope = -12 * ones(size(magfft,1),1); % -12 dB/octave
harmAmps =  zeros(size(magfft,1), numMaxHarmonics);
harmIdxs =  zeros(size(magfft,1), numMaxHarmonics);


maxPeakWin =  floor(4 * specsize / 2048); %4; if FFTsize is 4096, peakwidth is 8 bins
minPeakAmp = 10^(-90/20); % -100 dB
nHistAmpsCounter = 0;
for i=1:size(magfft,1)    
    if pitch(i)> 0 % 60Hz
        harmIdx = 1 + round([pitch(i): pitch(i):Fs/2]*freq2bin); % 0Hz correposnd to 1 for matlalb indexing
        harmIdx = harmIdx(1:min(length(harmIdx),numMaxHarmonics));
        idx = zeros(size(harmIdx));
        ampharm = zeros(size(harmIdx));
        
        % get peak(max) index for each harmonic                
        lbin = 1;
        rbin = 1;
        for k=1:length(harmIdx)-1,
            maxPeakWin2 = max(1, min(maxPeakWin, floor(0.5*(harmIdx(k+1) - harmIdx(k)))));
            lbin = max(rbin, harmIdx(k) - maxPeakWin2);
            rbin = min(harmIdx(k+1), harmIdx(k) + maxPeakWin2 -1);
            lbin = max(lbin,1);
            rbin = min(rbin,maxFreqBin);  
            % store default values
             idx(k) = harmIdx(k);
             ampharm(k) = magfft(i,harmIdx(k)); 
             maxReg = max(magfft(i,lbin:rbin));
             minReg = min(magfft(i,lbin:rbin));
            if( max(magfft(i,lbin:rbin)) > minPeakAmp) && ((rbin - lbin) > 4) && ((maxReg/minReg)  > 1.5)
                % TODO: find peaks is too intensive. Replace with faster
                % solution.
                %disp([i,k])
                [mm xx] = findpeaks(magfft(i,lbin:rbin));
                [msg id] = lastwarn;
                if size(id)>0,
                 warning('off',id); % set warning off for findpeaks which are annoying
                end
                
                if size(mm)> 0
                    [mm2 idxMax] = max(mm);                    
                    idx(k) = xx(idxMax) + lbin -1;
                    ampharm(k) = mm2;                    
                else                   
%                    disp('peaks not found');
                end  
            end % min peak amp
        end
        % smooth amplitude values in time for smooth harmonic envelope
        if i>1,
            %ampharm = 0.1 * ampharm + 0.9 * harmAmps(i-1,1:length(harmIdx));            
            nHistAmpsCounter = nHistAmpsCounter + 1; % update hisotry count
            ampsCounter = min(nHistAmpsCounter,10);
            ampharm = (1/ampsCounter) * ampharm + (1 - (1/ampsCounter)) * harmAmps(i-1,1:length(harmIdx));                        
        end
        
        % interpolate on the log scale
        harmIdxValid =  idx(idx>0);        % Keep only valid harmonic partials (non-zero indices) 
        ampharmValid = ampharm(idx>0);
        if size(harmIdxValid>0)
            %magDB = 20*log10(magfft(i,:));
            ampDB = 20*log10(ampharmValid); 
            itype = 'linear'; %'spline'
            if harmIdxValid(1)==1,            
                %f = interp1([harmIdxValid maxFreqBin],[magDB(harmIdxValid), min(magDB(harmIdxValid))],1:maxFreqBin,itype);
                  f = interp1([harmIdxValid maxFreqBin],[ampDB, min(ampDB)],1:maxFreqBin,itype);
            else
              %  f = interp1([1 harmIdxValid maxFreqBin],[magDB(harmIdxValid(1)), magDB(harmIdxValid), min(magDB(harmIdxValid))],1:maxFreqBin,itype);
                f = interp1([1 harmIdxValid maxFreqBin],[ampDB(1), ampDB, min(ampDB)],1:maxFreqBin,itype);
            end
            flin = 10.^(f/20);
            magfftEnv(i,1:maxFreqBin) = [flin.*(flin>0)];
        else
            magfftEnv(i,1:maxFreqBin) = magfft(i,1:maxFreqBin);
        end
    
        % estimate excitation slope
        excitSlope(i) = estExcitationSlope(magfftEnv(i,1:maxFreqBin), Fs);
        % force a fixed excitation slope above forceExSlopeFreq if
        % input envelope is higher.
        if (forceExSlopeFreq>0) 
            bin5K = min(floor(forceExSlopeFreq * freq2bin), size(magfftEnv,2)-30);
            meanSpec5K = mean(magfftEnv(i,bin5K-30:bin5K+30));
            const5K = (bin5K) ^(-excitSlope(i) / (20*log10(2)));
            slopeEnv = meanSpec5K * const5K * (bin5K:maxFreqBin).^(excitSlope(i) /(20*log10(2)));
            magfftEnv(i,bin5K:maxFreqBin) = min(magfftEnv(i,bin5K:maxFreqBin) , slopeEnv);
        end
        
        % store harmonic amplitudes
        harmAmps(i,1:length(harmIdx)) = ampharm;
        harmIdxs(i,1:length(harmIdx)) = idx;

     else  % --> now for unvoiced frames we compute as well the envelope
%     at equally-spaced bins 
          magfftEnv(i,1:maxFreqBin) =  zeros(1,maxFreqBin);
%          % copy at spaced intervals depending on the number of maximum
%          % harmonicx
%          magfftEnv(i,1:maxFreqBin) =  magfft(i,1:maxFreqBin);
%          harmIdx = 1:floor(maxFreqBin/maxNumHarmonics):maxFreqBin;
%         harmIdx = harmIdx(1:min(length(harmIdx),numMaxHarmonics));
%         ampharm = magfft(i,harmIdx);
%         harmAmps(i,1:length(harmIdx)) = ampharm;
%         harmIdxs(i,1:length(harmIdx)) = harmIdx;
        nHistAmpsCounter =  0 ; % reset the number of history frames for voiced amplitudes
     end
    
end

end

function  excitSlope = estExcitationSlope(magfft, Fs)
% init vars
excitSlope = -12; 
specsize = size(magfft,2);
freq2bin = 2*specsize/Fs;
% slope is estimated in the range 200Hz and 4kHz
freqs = 1 + floor(freq2bin * 200):floor(freq2bin * 4e3);
%mags = mean([20*log10(magfft(freqs)) ; 20*log10(magfft(freqs+3)) ; 20*log10(magfft(freqs-3))]);
p = polyfit( log10(freqs)./log10(2),20*log10(magfft(freqs)), 1);

excitSlope = p(1);
end
