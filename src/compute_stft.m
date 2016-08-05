function [magfft phasefft] = compute_stft(x, window, hopsize, fftsize, zerophase)
% function [magfft phasefft] = compute_stft(x, window, hopsize, fftsize)
% STFT analysis
% Jordi Janer, UPF 2010-2013.

% init vars
winsize = length(window);

doZeroPhase = zerophase;
windowNorm = sum(window)/2; % normaliza by window energy
% loop for all audio frames

x = [zeros(1,floor(winsize/2)) x zeros(1,floor(winsize/2))]; % add half window for first frames
lenX = length(x);

% init  matrices
magfft = [];%zeros(num_frames,winsize);
phasefft = [];%zeros(num_frames,winsize);
num_frames = size([0:hopsize:lenX-winsize],2);
magfft = zeros(num_frames,floor(fftsize/2) +1);
phasefft = zeros(num_frames,floor(fftsize/2) +1);


idxcount = 1;
for (i=0:hopsize:lenX-winsize),
    % center windowed signal and add zero padding
    win =  x(i+1 : i + (winsize) ).*window;
    centerPos = i + (winsize/2);

    % add zero-pad
    if doZeroPhase,
        win = [win(winsize/2+1:winsize) zeros(1,(fftsize - winsize)), win(1:winsize/2)]; % zero-phase
    else
        win = [win, zeros(1,fftsize - winsize)]; % nonzero phase
    end
    
    tf =  fft(win);
    tp = angle(tf);   
    
    %magfft = [magfft; abs(tf(1:floor(fftsize/2) +1)) / windowNorm]; % normalize spectrum
    %phasefft = [phasefft; tp(1:floor(fftsize/2) +1)];    
    magfft(idxcount,:) = [abs(tf(1:floor(fftsize/2) +1)) / windowNorm]; % normalize spectrum
    phasefft(idxcount,:) = [tp(1:floor(fftsize/2) +1)];

    idxcount = 1+idxcount;
end
end