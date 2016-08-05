function [magfft phasefft params] = get_stft(audiofile, inParams)
% 
% function [magfft phasefft] = get_stft(audiofilename, mono)
% It returns complex spectrogram (mono)
%
% Input: 
%   audio filename WAV
%   params: struct containing 
%       params.winsize, params.hopsize, params.fftsize (in samples)
% Output: 
%   magfft: magnitude spectrum
%   phasefft: phase spectrum

stftfilename = strrep(audiofile,'.wav',  '.stfts'); % HDF% files with STFT



% load STFT file if exists
if exist(stftfilename),
    disp(['Reading from: ', stftfilename]);
    
    if ~exist('inParams','var')
        params = {};
    end
    % using low-level Matlab funcitons to read HDF5 files
    % High-level functions do not work properly on Matlab 2009    
    hinfo =  hdf5info(stftfilename);
    fileid = H5F.open(stftfilename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
    dataID = H5D.open(fileid,'/stfts');    
    data = H5D.read(dataID, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');    
    paramID = H5D.open(fileid, '/parameters/sampleRate');%hinfo.GroupHierarchy.Groups.Datasets(8).Name);    
    params.fs = double(H5D.read(paramID, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT'));
    paramID = H5D.open(fileid, '/parameters/frameSize');%hinfo.GroupHierarchy.Groups.Datasets(8).Name);    
    params.winsize = double(H5D.read(paramID, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT'));            
    paramID = H5D.open(fileid, '/parameters/hopSize');%hinfo.GroupHierarchy.Groups.Datasets(8).Name);    
    params.hopsize = double(H5D.read(paramID, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT'));            
    paramID = H5D.open(fileid, '/parameters/fftSize');%hinfo.GroupHierarchy.Groups.Datasets(8).Name);    
    params.fftsize = double(H5D.read(paramID, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT'));            
    paramID = H5D.open(fileid, '/parameters/windowType');%hinfo.GroupHierarchy.Groups.Datasets(8).Name);    
    params.windowtype = (H5D.read(paramID, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT'))';            
    paramID = H5D.open(fileid, '/parameters/window');%hinfo.GroupHierarchy.Groups.Datasets(8).Name);    
    params.window = double(H5D.read(paramID, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT'))';            
    paramID = H5D.open(fileid, '/parameters/zeroPhase');
    params.zerophase = (H5D.read(paramID, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT'))';            
    H5F.close(fileid);
    
    magfft = abs(data.r + i * data.i)';
	phasefft = angle (data.r + i * data.i)';    
end 
if ~exist(stftfilename) || (exist('inParams.zerophase') && inParams.zerophase ~= params.zerophase),
    %disp(['Recomputing: ', stftfilename]);

    if ~exist('inParams','var') || ~isfield(inParams,'winsize')
        params.hopsize = 512;
        params.winsize = 4096;
        params.fftsize = 4096;
        params.zerophase = 1;
        params.windowtype = 'blackmanharris';
        sprintf ('Setting default STFT parameters: winsize %d, hopsize %d, wintype %s',  params.hopsize,  params.winsize,  params.windowtype)
    else
        params = inParams;
    end
    % read
    [x params.fs] = wavread(audiofile);        
    if (size(x,2) == 2), % stereo        
        x = 0.5 * [x(:,1) +  x(:,2) ];        
    end    
    
    window = hann(params.winsize)'; % default
    % create window Blackman-Harris 92dB
    if strcmp(params.windowtype,'blackmanharris')
        n = 0:params.winsize-1;
        BlackHar = 0.35875 - 0.48829*cos(2*pi*n/params.winsize) ...
            + 0.14128*cos(2*2*pi*n/params.winsize) - 0.01168*cos(3*2*pi*n/params.winsize);
        window = BlackHar;
    end
    if strcmp(params.windowtype,'hann')
        window = hann(params.winsize)';
    end
    [magfft phasefft] = compute_stft(x(:,1)', window, params.hopsize, params.fftsize, params.zerophase);
    params.window = window;
end


