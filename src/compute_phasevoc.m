function [syn_magfft, syn_phasefft, syn_pitch] = compute_phasevoc(magfft,phasefft, params, pitch, synPS, synTimeIdx, timbrePreserve)
% function [syn_magfft, syn_phasefft] = compute_phasevoc(magfft,phasefft,
% params, pitch, synPitch, synTimeIdx, timbrePreserve)
%
% It computes pitch-shift and time-stretch transformation for a monophonic sound. 
% Inputs:
%   - magfft: input magnitude spectrogram
%   - phasefft: input phase spectrogram
%   - params: structure with the FFT and algorithm configuration parameters
%   (hopsize, window_size)
%   - pitch: vector containign the pitch values (in Hz) per frame of the
%   input sound.
%   - synPS: vector containign the pitch shifting factor  (pitch_out / pitch_in) values per frame of the transformed  sound (output size).
%   - synTimeIdx: vector containing the time indices (in frames) of the
%   input sound used to generate the time-stertched output sound. (output size)
%   - timbrePreserve: flag to perserve the timbre of the input sound when
%   applying pitch-shifting.
%   
% outputs:

% ---------------------------------------------
% Init output spectrogram
syn_magfft = zeros(size(magfft)); 
syn_phasefft = zeros(size(phasefft));

if size(pitch,1) == 0,
    sprintf('Pitch data is empty.');
    return 
end

% ---------------------------------------------
% Configuration parameters
global useTimbrePreserve;
useTimbrePreserve = timbrePreserve; % for test sounds (e.g pure tone), timbre presenrvation will not work properly
global useFlatPhase;
useFlatPhase = false; %false; % uses flat phase for the harmonic peaks. It generates robotic artifacts
useInputPhase = true; % do not consider input phase information (or it does not exist).

% apply minipum pitch to process unvoiced frames with a fixed
pitch = pitch .* (pitch >= params.minPitchHz); 

% ---------------------------------------------
% Compute harmonic envelope, harmonic partials amplitudes and frequency bins indexes.
[magfftEnv excitSlope harmAmps harmIdxs] = harmonicEnv(magfft, pitch, params.fs, params.fs/2, params.numMaxHarmonics);
% compute Minimum phase values from harmonic partials 
minPhaseHarmonics = angle(mps(harmAmps)); % check number of harmonics. Default = 200

% compute output phase spectrum by interpolating/translating with original phase
inphasefft = zeros(size(magfft));
if useInputPhase
    inphasefft = phasefft;
end

% compute ouput spectrum from harmonic envelope
% single function for magnitude and phase
%--[syn_magfft, syn_phasefft, syn_pitch] = generateTransformedSpectra(magfft, magfftEnv, inphasefft, minPhaseHarmonics, pitch, harmIdxs, params.fs, params.hopsize, synPitch, synTimeIdx);
[syn_magfft, syn_phasefft, syn_pitch] = generateTransformedSpectra(magfft, magfftEnv, inphasefft, minPhaseHarmonics, pitch, harmIdxs, params, synPS, synTimeIdx);


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function generateTransformedSpectra()
%
% Generates the phase sepactrum from the pitch, minimum phase and original
% analysis spectrum
%--function [outmagfft, outphasefft, outpitch] = generateTransformedSpectra(magfft, magfftEnv, inphasefft, minPhaseHarmonics, pitch, harmIdxs, srate, hsize, synPitch, synTimeIdx)
function [outmagfft, outphasefft, outpitch] = generateTransformedSpectra(magfft, magfftEnv, inphasefft, minPhaseHarmonics, pitch, harmIdxs, params, synPS, synTimeIdx)
global useFlatPhase;
global useTimbrePreserve; 

% get params
srate = params.fs;
hsize = params.hopsize;

nFrAnal = size(inphasefft,1);
nFrSyn = size(synTimeIdx,1);

% init put buffers
outmagfft = zeros(nFrSyn,size(magfft,2));
outphasefft = zeros(nFrSyn,size(magfft,2));
outmpharms = zeros(size(minPhaseHarmonics));
outpitch = zeros(nFrSyn,1);

periodPos = 0; % init 
crossfadeFreq = 12e3;
specsize = size(inphasefft,2);
peakWidthBins = floor(4 * size(inphasefft,2) / 2049); % default value is 4 bins for a fftisze of 4096. 
freq2bin = 2*specsize/srate;

for j=1:nFrSyn, % loop for all synthesis frames
      
    % get analysis frame index
    i = round(synTimeIdx(j));        
    i = min(i, nFrAnal); 
    
    diffInterp  = synTimeIdx(j) - i;
    
    % compute the target minium phase for all harmonics with pitch offset
    mp = minPhaseHarmonics(i,:);  
    mp2 = zeros(size(mp));

    % Get analysis harmonic partials
    validIdx = find(harmIdxs(i,:));
    numValidHarmonics = 0;
    if size(validIdx)>0
        numValidHarmonics = validIdx(end);
    end
  %disp(num2str(100*j/nFrSyn))
    if pitch(i) > 0  && numValidHarmonics>0, %minPitchHz, 

            %analHarmIdxs = harmIdxs(i,1:numValidHarmonics);
           %inharmonicFactors = (analHarmIdxs * srate / (2*specsize))/pitch(i); % to deal with inharmonic patrials as in piano. It does not work properly, as a specific model in the peak selection is needed.
           
           % compute interpolated synthesis pitch value, usefule for long
           % time-stretch factors. Skip for first akd last frames
           % TODO: this needs to be checked as it produces strange modulations
            %--outpitch(j) = synPitch(i); % old
            outpitch(j) = pitch(i) * synPS(j); % new
            magInterp = magfft(i,:);
            magEnvInterp =  magfftEnv(i,:);
            phaseInterp = inphasefft(i,:);
            mpInterp = minPhaseHarmonics(i,:); 
            harmIdxsInterp = harmIdxs(i,1:numValidHarmonics);            
%--            synHarmIdxsInterp = harmIdxs(i,1:numValidHarmonics) * (synPitch(i)/ pitch(i));
%--            synHarmIdxCurFr = 1 + round([synPitch(i): synPitch(i):srate/2]*freq2bin); % 0Hz correposnd to 1 for matlalb indexing
            synHarmIdxsInterp = harmIdxs(i,1:numValidHarmonics) * (synPS(j)); % new
            synHarmIdxCurFr = 1 + round([pitch(i)*synPS(j):pitch(i)*synPS(j):srate/2]*freq2bin); % new -0Hz correposnd to 1 for matlalb indexing

            synHarmIdxCurFr = [synHarmIdxsInterp, synHarmIdxCurFr(numValidHarmonics+1:end)];
            useInterp = true; %true;%true; % true;
            if useInterp
                if diffInterp <0 && (i > 1)
                        % synthesis index init
%--                        synHarmIdxLastFr = 1 + round([synPitch(i-1): synPitch(i-1):srate/2]*freq2bin); % 0Hz correposnd to 1 for matlalb indexing                        
%--                        synHarmIdxLastFr = [harmIdxs(i-1,1:numValidHarmonics) * (synPitch(i-1)/pitch(i-1)), synHarmIdxLastFr(numValidHarmonics+1:end)];     
                        synHarmIdxLastFr = 1 + round([pitch(i-1)*synPS(j):pitch(i-1)*synPS(j):srate/2]*freq2bin); % 0Hz correposnd to 1 for matlalb indexing                        
                        synHarmIdxLastFr = [harmIdxs(i-1,1:numValidHarmonics) * (synPS(j)), synHarmIdxLastFr(numValidHarmonics+1:end)];     
                        
                        synHarmIdxCurFr = synHarmIdxCurFr(1:min(size(synHarmIdxCurFr,2),size(synHarmIdxLastFr,2))); 
                        synHarmIdxLastFr = synHarmIdxLastFr(1:size(synHarmIdxCurFr,2));
                    if  pitch(i-1)>0,
    %--                    outpitch(j) = (-diffInterp) * synPitch(i-1) + synPitch(i) *(1+ diffInterp);
                        outpitch(j) = (-diffInterp) * pitch(i-1)*synPS(j) + pitch(i) * synPS(j) *(1+ diffInterp);
                        magInterp = (-diffInterp) * magfft(i-1,:) + magfft(i,:) *(1+ diffInterp);
                        magEnvInterp = (-diffInterp) * magfftEnv(i-1,:) + magfftEnv(i,:) *(1+ diffInterp);
                        phaseInterp = (-diffInterp) * inphasefft(i-1,:) + inphasefft(i,:) *(1+ diffInterp);
                        mpInterp = (-diffInterp) * minPhaseHarmonics(i-1,:) + minPhaseHarmonics(i,:) *(1+ diffInterp);
                        harmIdxsInterp = (-diffInterp) * harmIdxs(i-1,1:numValidHarmonics) + harmIdxs(i,1:numValidHarmonics) *(1+ diffInterp);                        
                        %synHarmIdxsInterp = (-diffInterp) * harmIdxs(i-1,1:numValidHarmonics) * (synPitch(i-1)/pitch(i-1)) + harmIdxs(i,1:numValidHarmonics) * (synPitch(i)/pitch(i))*(1+ diffInterp);                        
                        synHarmIdxsInterp = (-diffInterp) * synHarmIdxLastFr + synHarmIdxCurFr *(1+ diffInterp);                        
                    end
                end
                if  diffInterp >=0 && (i < nFrAnal)
%--                    synHarmIdxNextFr = 1 + round([synPitch(i+1): synPitch(i+1):srate/2]*freq2bin); % 0Hz correposnd to 1 for matlalb indexing
%--                    synHarmIdxNextFr = [harmIdxs(i+1,1:numValidHarmonics) * (synPitch(i+1)/pitch(i+1)), synHarmIdxNextFr(numValidHarmonics+1:end)];
                    synHarmIdxNextFr = 1 + round([pitch(i+1)*synPS(j): pitch(i+1)*synPS(j):srate/2]*freq2bin); % 0Hz correposnd to 1 for matlalb indexing
                    synHarmIdxNextFr = [harmIdxs(i+1,1:numValidHarmonics) * (synPS(j)), synHarmIdxNextFr(numValidHarmonics+1:end)];
                    
                    synHarmIdxCurFr = synHarmIdxCurFr(1:min(size(synHarmIdxCurFr,2),size(synHarmIdxNextFr,2))); 
                    synHarmIdxNextFr = synHarmIdxNextFr(1:size(synHarmIdxCurFr,2));
                    if  pitch(i+1)>0,
                        %-- outpitch(j) = (1-diffInterp) * synPitch(i) + synPitch(i+1) * (diffInterp);
                        outpitch(j) = (1-diffInterp) * pitch(i)*synPS(j) + pitch(i+1)*synPS(j) * (diffInterp);
                        magInterp = (1-diffInterp) * magfft(i,:) + magfft(i+1,:) * (diffInterp);
                        magEnvInterp = (1-diffInterp) * magfftEnv(i,:) + magfftEnv(i+1,:) * (diffInterp);                        
                        phaseInterp = (1-diffInterp) * inphasefft(i,:) + inphasefft(i+1,:) * (diffInterp);                        
                        mpInterp = (1-diffInterp) * minPhaseHarmonics(i,:) + minPhaseHarmonics(i+1,:) * (diffInterp);
                        harmIdxsInterp = (1-diffInterp) * harmIdxs(i,1:numValidHarmonics) + harmIdxs(i+1,1:numValidHarmonics) * (diffInterp);
                        %synHarmIdxsInterp = (1-diffInterp) * harmIdxs(i,1:numValidHarmonics) * (synPitch(i)/pitch(i)) + harmIdxs(i+1,1:numValidHarmonics) * (synPitch(i+1)/pitch(i+1)) * (diffInterp);
                        synHarmIdxsInterp = (1-diffInterp) * synHarmIdxCurFr + synHarmIdxNextFr * (diffInterp);
                    end
                end
            end % use interpolation
            
            
           % Compute normalized period position (0,1) according to current
           % pitch value and pitch-shift factor. Take into account if
           % time-expansio to average pitchvalues or take
           % interpolated synthesis pitch values.            
           PSfactor =  outpitch(j) / pitch(i);   % synPitch(i)        
           period = 1 / (pitch(i)*PSfactor);  

           % interpolated output pitch values
           if j>1 && outpitch(j-1)>0,
               period = 1 / (mean(outpitch(j-1:j)));      
                PSfactor =  mean(outpitch(j-1:j)) / pitch(i); 
           end
           
           normalizedPos = mod(periodPos + ((hsize/srate) / period) ,1);
           periodPos = normalizedPos;           
           
           % Round the interpolated frequnecy bin indexes for synthesis harmonic partials           
           synHarmIdxs = round(synHarmIdxsInterp); % new
           synHarmIdxs = min(synHarmIdxs, specsize);
           diffHarmIdxs = synHarmIdxsInterp - synHarmIdxs; 
           analHarmIdxs = round(harmIdxsInterp); % new use interpolation indices
           analHarmIdxs = min(analHarmIdxs, specsize);
           diffAnalHarmIdxs = harmIdxsInterp - analHarmIdxs;

           
           % Initialize lower bins. Set random phase for the first  bins until half the first harmonic
           l = 1;
           r = max(1,floor((analHarmIdxs(1)/2))); % begin bin.
           l2 = 1;
           r2 =  max(1,floor((synHarmIdxs(1)/2))); % begin bin.
           
           % Fill first bins of the spectrum until first harmonic partial.
           outmagfft(j,l2:r2) =  10^(-90/20)*magfft(i,l2:r2); % copy attenuated frist input bins
           outphasefft(j,l2:r2) = 2*pi*(-0.5 + rand(size(inphasefft(i,l2:r2))));            
           l2 = floor((synHarmIdxs(1)/2)); 
           r2 = l2;
           peakWidth = floor(min(2*specsize/((PSfactor*period*srate)), 2*specsize/((period*srate))));
           
           % Compute the output phase for all synthesis partials
           for k = 1:size(synHarmIdxs,2),
               
               % find closest analysis partial index for each synthesis partial
               % n: index of the analysis partial bins
               % k: index of the synthesis partial bins
               tmp = abs(analHarmIdxs-synHarmIdxs(k));
               [midx closestIdx] = min(tmp); %index of closest value                 
               n = closestIdx;   
               
               % find second closest peak to interpolate
               n_2nd = min(size(harmIdxsInterp,2),(n+1));
               if n>1 && (n <size(harmIdxsInterp,2))
                   if abs(harmIdxsInterp(n-1) - analHarmIdxs(n)) > abs(harmIdxsInterp(n+1) - analHarmIdxs(n))
                       n_2nd = min(size(harmIdxsInterp,2),(n+1));
                   else
                       n_2nd = (n-1);
                   end
               end
               % distance to harmonic analysis indexes for interpolation  
               dist_n1st = abs(harmIdxsInterp(n)-synHarmIdxs(k));
               dist_n2nd = abs(harmIdxsInterp(n_2nd)-synHarmIdxs(k));
               dist_nTotal = dist_n1st + dist_n2nd;
               
               if synHarmIdxs(k) == 0 || analHarmIdxs(n) == 0,
                   continue;
               end
                           
               % find leftindexes 
               l =  max(1,min(specsize,analHarmIdxs(n) - floor(peakWidth/2)));
               l2 = max(1,min(specsize,synHarmIdxs(k) - floor(peakWidth/2))); 
               l_2nd =  max(1,min(specsize,analHarmIdxs(n_2nd) - floor(peakWidth/2)));
               peakWidth = min(peakWidth,specsize - l2); % avoid out of bounds index
               peakWidth = min(peakWidth,specsize - l); 
                
       
        
               % Estimate new phase
               % compute minimum-phase with partial continuation
               %mp2(k) = mp(n) + 2*pi*(k)*periodPos;
               mp2(k) = mpInterp(n) + 2*pi*(k)*periodPos;
               %mp2(k) = mp(n) + 2*pi*(inharmonicFactors(k))*periodPos; % taking into account non-harmonic relations of patrials (e.g. in paino)
               
               % Shift the phase of all bins for that partial
               inPhase = phaseInterp(analHarmIdxs(n)); %  phase value at analysis peak center
               
               
                % Replace peak magnitude (initialized peak copy)
                outmagfft(j,l2:l2+peakWidth) = magInterp(l:l+peakWidth); % using interpolated spectrum                
                
                % --------------------------------
                % Interpolate analysis peak according to analysis partial bin indices
                magInterpPeak = magInterp(l:l+peakWidth);
                phaseInterpPeak = phaseInterp(l:l+peakWidth);
                % changing k to n for indexing diffAnalHarmIdxs / 21-5-14
                if diffAnalHarmIdxs(n) < 0 && (l>1)
                    magInterpPeak = (-diffAnalHarmIdxs(n)) * magInterp(l-1:l-1+peakWidth) + (1 + diffAnalHarmIdxs(n)) * magInterp(l:l+peakWidth);
                    magInterpPeak2nd = (-diffAnalHarmIdxs(n_2nd)) * magInterp(l_2nd-1:l_2nd-1+peakWidth) + (1 + diffAnalHarmIdxs(n)) * magInterp(l_2nd:l_2nd+peakWidth);
                    phaseInterpPeak = (-diffAnalHarmIdxs(n)) * phaseInterp(l-1:l-1+peakWidth) + (1 + diffAnalHarmIdxs(n)) * phaseInterp(l:l+peakWidth);
                end
                if diffAnalHarmIdxs(n) >= 0 && (l>1)
                    magInterpPeak = (1-diffAnalHarmIdxs(n)) * magInterp(l:l+peakWidth) + (diffAnalHarmIdxs(n)) * magInterp(l+1:l+1+peakWidth);
                    magInterpPeak2nd = (1-diffAnalHarmIdxs(n_2nd)) * magInterp(l_2nd:l_2nd+peakWidth) + (diffAnalHarmIdxs(n_2nd)) * magInterp(l_2nd+1:l_2nd+1+peakWidth);
                    phaseInterpPeak = (1-diffAnalHarmIdxs(n)) * phaseInterp(l:l+peakWidth) + (diffAnalHarmIdxs(n)) * phaseInterp(l+1:l+1+peakWidth);
                end
                
                % interpolation to avoid jumps in the harmonic analysius
                % index on consecutive frames.
                if (dist_nTotal>0)
                    magInterpPeak = (dist_n2nd/dist_nTotal) * magInterpPeak + (dist_n1st/dist_nTotal)*magInterpPeak2nd;
                end
                % --------------------------------
                % output phase per partial               
                if useFlatPhase,
                    outphasefft(j,l2:l2+peakWidth) = -(inPhase - mp2(k));
                else
%                     outphasefft(j,l2:l2+peakWidth) = phaseInterp(l:l+peakWidth) - (inPhase - mp2(k));
%                     % Interpolate copied synthesis peak  according to
                    % synthesis partial bin indices  
                    if diffHarmIdxs(k) <0
                        outphasefft(j,l2:l2+peakWidth)  = ((1+diffHarmIdxs(k)) * phaseInterpPeak + ...
                            (-diffHarmIdxs(k)) * [phaseInterpPeak(2:2+peakWidth-1) 0] ) - (inPhase - mp2(k)); % using interpolated spectrum
                    end
                    if (diffHarmIdxs(k) >= 0)  && (l > 1)
                        outphasefft(j,l2:l2+peakWidth)  = ((1-diffHarmIdxs(k)) * phaseInterpPeak + ...
                            (diffHarmIdxs(k)) * [0 phaseInterpPeak(1:1+peakWidth-1)]) - (inPhase - mp2(k)) ; % using interpolated spectrum
                    end 
                    
                end
               
                % --------------------------------                           
                % output magnitude per partial
                % Interpolate copied synthesis peak  according to synthesis partial bin indices
                if diffHarmIdxs(k) <0
                    outmagfft(j,l2:l2+peakWidth) = (1+diffHarmIdxs(k)) * magInterpPeak + ...
                        (-diffHarmIdxs(k)) * [magInterpPeak(2:2+peakWidth-1) 0] ; % using interpolated spectrum
                end
                if (diffHarmIdxs(k) >= 0)  && (l > 1)
                    outmagfft(j,l2:l2+peakWidth) = (1-diffHarmIdxs(k)) * magInterpPeak + ...
                        (diffHarmIdxs(k)) * [0 magInterpPeak(1:1+peakWidth-1)] ; % using interpolated spectrum
                end
                
                % --------------------------------
                % Post processing of output magnitude and phases 
                % Fill the gaps  between synthesis partials if necessary              
                if k>1 &&  (synHarmIdxs(k-1) > 0)
                    r2prev = min(synHarmIdxs(k-1) + floor(peakWidth/2),specsize);
                    leftPeakDiff =  20*log10( outmagfft(j,synHarmIdxs(k-1)) / (outmagfft(j,r2prev)));
                    rightPeakDiff = 20*log10( outmagfft(j,synHarmIdxs(k)) / (outmagfft(j,l2)));
                    diffDB = 20;
                    % interpolate from partial bounds if left and right partials have a peak
                    % shape to avoid discontinuities in the spectrum producing time aliasing
                    if l2>r2prev && leftPeakDiff > diffDB && rightPeakDiff > diffDB, 
                        outmagfft(j,r2prev:l2) = interp1([r2prev, l2],[outmagfft(j,r2prev), outmagfft(j,l2)],[r2prev:l2],'linear');
                    else
                        outmagfft(j,r2prev:l2) = 10^(-90/20)*outmagfft(j,synHarmIdxs(k));
                    end
                end     
                                              
                % add random noise modifying the phase for upper harmonics (> crossFadeFreq)             
                if (k*pitch(i)*PSfactor > crossfadeFreq),
                    outphasefft(j,l2:l2+peakWidth) =  (inPhase - mp2(k)); %  force flat phase for all bins in the harmonic                                         
                    phaseInterpolate = hann(peakWidth+1)';
                    outphasefft(j,l2:l2+peakWidth) = outphasefft(j,l2:l2+peakWidth) .*  phaseInterpolate + (1- phaseInterpolate) .* (2*pi*(-0.5 + rand(1,peakWidth+1)));                       
                end
                % --------------------------------
                
                        
                
           end    % loop for all harmonics
           
           
        % Preserve timbre envelope
        if (useTimbrePreserve && size(synHarmIdxs,2) > 3)
            % compute filter to preserve input timbre. Note that for test
            % signals (e.g. pure tones), the harmonic envelope might be
            % invalid.          
            %envFilter = computeTimbreFilter(magfftEnv(i,:), outmagfft(j,:), pitch(i)*PSfactor, srate);        
            envFilter = computeTimbreFilter(magEnvInterp, outmagfft(j,:), pitch(i)*PSfactor, params);              
            outmagfft(j,:) = outmagfft(j,:) .* envFilter;       
        end
        
        % Apply gain factor to keep original frame energy 
        gainFactor = sqrt(sum(magInterp.^2)/ sum(outmagfft(j,:).^2)); %   ENERGY 
       %  gainFactor = sum(magInterp)/ sum(outmagfft(j,:)); % MAG 
        outmagfft(j,:) = (gainFactor) * outmagfft(j,:);                      
        
        % End process for pitched frames
        % ----------------------------
    else
        % Process for unpitched frames (pitch < minPitchHz)
        % Magnitude and Phase spectra are copied from analysis frame.               
        periodPos = 0;
        TS_factor = 1;
        if j > 1
            TS_factor = (1 / (synTimeIdx(j) - synTimeIdx(j-1)));
        end
        if  TS_factor <= 1,
            outphasefft(j,:) = inphasefft(i,:); % use input phase if TS <= 1
        else
            outphasefft(j,:) =  rand(size(inphasefft(i,:)))*2*pi - pi;
        end
        
        % Magnitude interpolation for unvoiced frames
        magInterp = magfft(i,:);
        if diffInterp <0 && (i > 1)
            magInterp = (-diffInterp) * magfft(i-1,:) + magfft(i,:) *(1+ diffInterp);
        end
        if  diffInterp >=0 && (i < nFrAnal)            
            magInterp = (1-diffInterp) * magfft(i,:) + magfft(i+1,:) * (diffInterp);            
        end        
        outmagfft(j,:) = magInterp;           
    end      
    % round  phase in range -pi, pi    
    outphasefft(j,:) =  roundphase(outphasefft(j,:)) ;
end


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From JOS matlab code : https://ccrma.stanford.edu/~jos
function [sm] = mps(magfft)
% [sm] = mps(s)
% create minimum-phase spectrum sm from complex spectrum s
sm = zeros(size(magfft));
for i=1:size(magfft,1)
    s = magfft(i,:); % frame
    s = [s s(end-1:-1:2)];
    smtemp = exp( fft( fold( ifft( log( clipdb(s,-100) )))));
    sm(i,:) = smtemp(1:size(magfft,2));
end
end


function [clipped] = clipdb(s,cutoff)
% [clipped] = clipdb(s,cutoff)
% Clip magnitude of s at its maximum + cutoff in dB.
% Example: clip(s,-100) makes sure the minimum magnitude
% of s is not more than 100dB below its maximum magnitude.
% If s is zero, nothing is done.

clipped = s;
as = abs(s);
mas = max(as(:));
if mas==0, return; end
if cutoff >= 0, return; end
thresh = mas*10^(cutoff/20); % db to linear
toosmall = find(as < thresh);
clipped = s;
clipped(toosmall) = thresh;
end

function [rw] = fold(r)
% [rw] = fold(r)
% Fold left wing of vector in "FFT buffer format"
% onto right wing
% J.O. Smith, 1982-2002

[m,n] = size(r);
if m*n ~= m+n-1
    error('fold.m: input must be a vector');
end
flipped = 0;
if (m > n)
    n = m;
    r = r.';
    flipped = 1;
end
if n < 3, rw = r; return;
elseif mod(n,2)==1,
    nt = (n+1)/2;
    rw = [ r(1), r(2:nt) + conj(r(n:-1:nt+1)), ...
        0*ones(1,n-nt) ];
else
    nt = n/2;
    rf = [r(2:nt),0];
    rf = rf + conj(r(n:-1:nt+1));
    rw = [ r(1) , rf , 0*ones(1,n-nt-1) ];
end;

if flipped
    rw = rw.';
end
end
% 



% function to wrap phase
function rphase = roundphase(ph)
rphase = mod(ph,2*pi);
for i=1:length(ph)
    if rphase(i)< -pi
        rphase(i) = rphase(i) + 2*pi;
    else if rphase(i)> pi
            rphase(i) = rphase(i) - 2*pi;
        end
    end
end
end
