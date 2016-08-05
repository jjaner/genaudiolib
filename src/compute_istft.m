
function [yout2] = compute_istft(outmagfft, outphasefft, window, hopsize, fftsize, overlapWinSize, type, zeroPhase)
%function [yout2] = compute_istft(outmagfft, outphasefft, window, hopsize, fftsize, zp, overlapWinSize, type, zeroPhase)
% inverse STFT
% Jordi Janer, UPF 2010.
% Params: 
%   - outmagfft, outphasefft: input spectrogram in magnitude and phase
%   - window: analysis window (vector)
%   - hopsize: analysis hopsize
%   - fftsize: analysis fftsize
%   - overlapSize: only used in triangular overlap-add
%   - type: 'tukey' or 'triangular'
%   - zeroPhase: spectrum has been cmputed using zero-phase


num_frames = size(outmagfft,1);
winsize = length(window);
yout = zeros(1,num_frames * hopsize + winsize * 2);
yout2= yout;


doZeroPhase = zeroPhase; %true;
    
% Select overlap window
    if winsize ~= fftsize
        disp('Error: Zero padding not yet implemented. Please use the option zero-padding = 0' )
    end

    yout2= yout;
    % crete overlap window
    if strcmp(type,'triangular'),
        overlapWin = triang(overlapWinSize+1)';
    end
    if strcmp(type,'tukey')
        overlapWin = tukeywin(overlapWinSize+1,1)';        
    end
    overlapWin = overlapWin(1:overlapWinSize);
    
    tk = [zeros(1,(winsize-overlapWinSize)/2),overlapWin,zeros(1,(winsize-overlapWinSize)/2)]; % original       
    
    overlapFactor = zeros(1,hopsize);
    for (i=1:hopsize:winsize-hopsize+1)        
           overlapFactor = overlapFactor +  window(i:i+hopsize-1).* tk(i:i+hopsize-1);
    end
    windowNorm = sum(window)/2;
    %for (i=winsize/(2*hopsize)+1: num_frames + winsize/(2*hopsize)),
    for (i=winsize/(2*hopsize)+1: num_frames + winsize/(hopsize)),
        % for (i=1:num_frames),
        
        j = i - winsize/(2*hopsize); % frame idx   
        if j >num_frames
            tf = zeros(1,fftsize);
        else
        tf = [windowNorm * outmagfft(j,1:fftsize/2 +1).*exp(1i*outphasefft(j,1:fftsize/2 +1)) , windowNorm * outmagfft(j,fftsize/2:-1:2).*exp(-1i*outphasefft(j,fftsize/2:-1:2))]; % denormalize spectrum (multiply by winsize)
        end
        
        ym = real(ifft(tf));
        ymz = ym;
        if doZeroPhase,
            ymz(1:winsize/2) = ym(winsize/2+1:winsize);
            ymz(winsize/2+1:winsize) =  ym(1:winsize/2);
        end
        ym = ymz(1:winsize) .* tk;        
        
        outindex =  [((i-1)*hopsize - (winsize/2) +1: (i-1)*hopsize + (winsize/2))];
        yout(outindex) = yout(outindex) + ym;
        outindex2 = [outindex(1):outindex(1) + hopsize-1];        
        yout2(outindex2) = yout(outindex2) ./ overlapFactor;                 
    end
    
%end,

% 
% if strcmp(type,'triangular'),
%     triang = [zeros(1,winsize/2 - (overlapWinSize/2)) , 0:1/(overlapWinSize/2):1, ...
%         1-(1/(overlapWinSize/2)) : -1/(overlapWinSize/2):1/(overlapWinSize/2), zeros(1,winsize/2 - (overlapWinSize/2))];
%     for (i=winsize/(2*hopsize)+1:num_frames-1),
%         j = i - winsize/(2*hopsize); % frame idx
%         tf = [outmagfft(j,1:fftsize/2 +1).*exp(1i*outphasefft(j,1:fftsize/2 +1)) , outmagfft(j,fftsize/2:-1:2).*exp(-1i*outphasefft(j,fftsize/2:-1:2))];
%         ym = real(ifft(tf));
%         ymz = ym;
%         if doZeroPhase,
%             ymz(1:winsize/2) = ym(winsize/2+1:winsize);
%             ymz(winsize/2+1:winsize) =  ym(1:winsize/2);
%         end
%         ym = ymz(1:winsize) ./ window;
%         
%         outindex = ((i-1)*hopsize - (winsize/2) +1: (i-1)*hopsize + (winsize/2));
%         yout2(outindex) = yout(outindex) + ym.*triang;
%     end
% end, % overlap window type

% adjust initial offset of half window
fistSample = (winsize/2)+1;
lastSample = (winsize/2) + num_frames*hopsize;
%yout2 = yout2((winsize/2)+1:end);
yout2 = yout2(fistSample:lastSample);
end