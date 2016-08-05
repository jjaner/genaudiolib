function apl_masking(audiofile, out_folder)
% function [W, H, outfile] = apl_masking(audiofile, nmf_model, param)
% SMC Audio Processing Lab 2015 - STFT-based SOFT MASKING SKELETON
%
% inputs: 
% - audiofile: input audio filename in WAV mono format.
% - out_folder: output filder where output separated files are copied
% 
% Jordi Janer, 2015, MTG-UPF
% 

%-----------------------------------------------------------
% Configuration parameters

% set initialization parametrers
param = apl_init();

% General: input file %
inputFileName = audiofile;


%-------------------------------------------------------------
% TIME-FREQUENCY ANALYSIS
%-------------------------------------------------------------
% applying stft to signals and obtaining an amplitude spectrogram.

[Y Yphase param] = get_stft(inputFileName, param);
nFrames = size(Y,1);

% initialize output soft masks for target and other to one (by-pass)
maskTarget = ones(size(Y'));
maskOther = ones(size(Y'));



%-------------------------------------------------------------
% CREATE SOFT MASK
%-------------------------------------------------------------
% add your code here to compute "target" and "other" soft masks for the
% separation
% maskTarget =  ...
% maskOther = ...

                                           
%-------------------------------------------------------------
% SIGNAL SEPARATION - SOFT MASKING AND RESYNTHESIS
%-------------------------------------------------------------

% Magnitude
outmagTarget = (maskTarget .* Y')';
outmagOther = (maskOther .* Y')';
% Phase
outphaseTarget = Yphase;
outphaseOther = Yphase;

% write output files
[~, audio_filename] = fileparts(audiofile);
output_files{1} = fullfile(out_folder,[audio_filename,'_target.wav']);
output_files{2} = fullfile(out_folder,[audio_filename,'_other.wav']);

write_stft(outmagTarget, outphaseTarget, output_files{1}, param)
write_stft(outmagOther, outphaseOther, output_files{2}, param)

                                           
                                           
         