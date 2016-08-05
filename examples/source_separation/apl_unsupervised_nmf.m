function apl_unsupervised_nmf(audiofile, out_matfile, out_folder)
% function [W, H, outfile] = apl_unsupervised_nmf(audiofile, nmf_model, param)
% SMC Audio Processing Lab 2015 - DECOMPOSITION STEP
%
% Descrition: 
%   This script estimates NMF basis learnt in an unsupervised fashion. It has the option to add additional constraints
%   on the activations matrix H0 in the context of informed source separation (e.g. MIDI is available).
%   The resulting matrix W can be stored as a .mat file, to be loaded and
%   reused in the apl_supervised_NMF() function to performed further separation.
%
% inputs: 
% - audiofile: input audio filename in WAV mono format.
% - out_matfile: output MAT file containing the W target matrix
% - out_folder: output filder where output audio files
%
%  
% Jordi Janer, 2015, MTG-UPF
% 

% set initialization parametrers
% Look at it to change parameters such as hopsize, framesize,etc.
param = apl_init();

% include paths
addpath('./3rdparty/genaudiolib');
addpath('./3rdparty/nmflib');
addpath('./3rdparty/midi_lib');


%-----------------------------------------------------------
% Configuration parameters

% General: input file %
inputFileName = audiofile;

%-------------------------------------------------------------
% TIME-FREQUENCY ANALYSIS
%-------------------------------------------------------------
% applying stft to signals and obtaining an amplitude spectrogram.

[Y, Yphase, param] = get_stft(inputFileName, param);
nFrames = size(Y,1);

% initialize output soft masks for target and other to one (by-pass)
maskTarget = ones(size(Y'));
maskOther = ones(size(Y'));



%-------------------------------------------------------------
% SIGNAL DECOMPOSITION - NMF
%-------------------------------------------------------------



% Initialize W_target with random values
n_basis = 50; % 50 basis for target component. (default: 50)
disp('TODO: check sizes for piano model to match thge unspuervsixed case')
W = rand(size(Y,1), n_basis); % matrix size is specsize x num_basis


%-------------------------------------------------------------
% Informed Source Separation (optional)
% Additional constrains for NMF estimation step

% H0 constraints for activations of target and other basis (e.g. for score-informed case)
H0mask = ones(size(W,2),nFrames); % set to 0 activations for W(target) that will be inactive. (default: 1)
% --> ADD YOUR CODE HERE TO INITIALIZE H0mask

% generate random initialization for target and other activation gains
nBasisConstrain = size(W,2);
H0init = [rand(size(H0mask)) .* H0mask];


%-------------------------------------------------------------
% Factorization

fprintf('*** Performing factorization ***\n');
[Wout,Hout,errs,vout] = nmf_kl_con(Y',nBasisConstrain, 'niter', param.niter, 'win', param.width,'W0',W,'H0',H0init,'norm_w',param.norm_w, 'n_other',param.n_other);
numBasisOut = size(Wout,2);

% Save W matrix to be used in a supervsied NMF 
save (out_matfile,'Wout'); % 




                                           
%-------------------------------------------------------------
% SIGNAL SEPARATION - SOFT MASKING AND RESYNTHESIS (optional)
%-------------------------------------------------------------


% define which basis are target and other to perform separation
targetBasisIdx = [1:nBasisConstraint]; % (default: 1:nBasisConstraint)
otherBasisIdx = [];
% --> ADD YOUR CODE HERE TO INITIALIZE TARGET AND OTHER BASIS INDEXES


% target basis/activations are Wt and Ht
Wt = Wout(:,targetBasisIdx,:);
Ht = Hout(targetBasisIdx,:);
% other basis/activations initilized to zero when not used
Wo = Wout(:,otherBasisIdx,:);
Ho = Hout(otherBasisIdx,:);

% reconstruct estimated target and other spectrogram
% convolutive NMF
target =  rec_cnmf(Wt,Ht,param.eps); % Wt * Ht;
other = rec_cnmf(Wo,Ho,param.eps); %Wo * Ho;

% compute filter to TF input signal (Wiener filter)
maskTarget = double(target) ./ (double(target) + double(other) + eps); % avoid division by zero
maskOther = double(other) ./ (double(target) + double(other) + eps);



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

                                           
                                           
         