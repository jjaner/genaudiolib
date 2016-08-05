function write_stft(magfft, phasefft, outfile, params)

% create output filder if needed
[path, name, ext ] = fileparts(outfile);
if ~exist(path,'dir')
    mkdir(path);
end

if ~isfield(params,'overlapWinSize')
    params.overlapWinSize = params.winsize;
end

% create analysis window 
window = params.window; % hann(params.winsize)'; % default
if isempty(window),
    if strcmp(params.windowType,'blackmanharris')
        %Blackman-Harris 92dB
        n = 0:params.winsize-1;
        BlackHar = 0.35875 - 0.48829*cos(2*pi*n/params.winsize) ...
            + 0.14128*cos(2*2*pi*n/params.winsize) - 0.01168*cos(3*2*pi*n/params.winsize);
        window =  BlackHar;
    end
    if strcmp(params.windowType,'hann')
        % hann window
        window = hann(params.winsize)';
    end
end
[yout] = compute_istft(magfft, phasefft, window, params.hopsize, params.fftsize, params.overlapWinSize, 'tukey', params.zerophase);
%[yout] = compute_istft(magfft, phasefft, window, params.hopsize, params.fftsize, params.overlapWinSize, 'triangular', params.zeroPhase);
wavwrite([yout'],params.fs,16,outfile);
