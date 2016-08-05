function concatenateBlocks(infolder, outfile, sourceName, hopsize, blocksize, overlapsize)
% function concatenateBlocks(infolder, outfile, sourceName, hopsize, blocksize, overlapsize)
% 
% Concatenation of audio files processed block-wise by source separation
% algorithms. Block audio files (WAV) need to be generated with overlap to be
% later concatenated. Block audio files should have this filename:
% SourcName_block_<BlockNumber>.wav
% 
% Inputs:
% - infolder: input  folder with the individual files. 
%             Input filename shall have the following format:
%             SourcName_block_<BlockNumber>.wav . BlockNumber starting at
%             0.
%             Input  audio files need to be in WAV format (mono, stereo).
% - outfile: output  filename
% - sourceName: instrument (source) name 
% - hopsize: hopsize in secs
% - blocksize: number of seconds contained in one block audio file
% - overlapsize: number of secs used to overlap
%
% 
% BLOCK CONCATENATION - MiesPV
% MTG-UPF 2014, Jordi Janer


blockfiles = dir(fullfile(infolder, [sourceName '_*.wav']));
nBlocks = size(blockfiles,1);

% get Fs and number of channels
[in Fs]= wavread(fullfile(infolder, blockfiles(1).name));
overlapsize=ceil(overlapsize*Fs);
hopsize=ceil(hopsize*Fs);

nCh = size(in,2);

x_out = zeros(nBlocks * blocksize * Fs, nCh);

% crossfade ramp
ramp = [0:1/(overlapsize-1): 1]';
if nCh ==2 ,
    ramp = [ramp , ramp];
end


for i=0:nBlocks-1,

    % find block filename per order
    blockidx = -1;
    iter=0;
    while blockidx ~= i,
        iter=iter+1;
        res = sscanf(blockfiles(iter).name,[sourceName,'_block_%d.wav']);
            blockidx = res(1);     
    end    
   
   infile = fullfile(infolder, blockfiles(iter).name);
   in = wavread(infile);    
   size_in = size(in,1);
       
   % merge crosssfade
   beginTime = blockidx * blocksize; % in seconds
   nbeg = 1 + floor(beginTime*Fs);   % In samples   
   %disp([nbeg, size_in])
   if i == 0 % 1
       x_out (1:size_in,:) =  in(1:size_in,:);
   else
       % overlap 
       overlapsize = min(overlapsize, size_in - (overlapsize -1));
       x_out (nbeg:nbeg+overlapsize-1,:) = (1- ramp(1:overlapsize)).* x_out(nbeg:nbeg+overlapsize-1  ,:) + ramp(1:overlapsize) .* in(1:overlapsize,:);
       % rest
       x_out (nbeg+overlapsize:nbeg+size_in-1,:) =  in(overlapsize+1:size_in,:);
   end
end

doNormalize = false;
if doNormalize,
    x_out=x_out./(1.1*max(abs(x_out)));
end
wavwrite (x_out,Fs, 16, outfile);