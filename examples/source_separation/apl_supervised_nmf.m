function apl_supervised_nmf(audiofile, out_folder)
% function [W, H, outfile] = apl_supervised_nmf(audiofile, nmf_model, param)
% SMC Audio Processing Lab 2015 - SUPERVISED SEPARATION STEP
%
% Descrition: 
%   This script performs source separation of an input audio file in two
%   components: "target" and "other". Target component is separated in a
%   supervised-NMF fashion for which a matrix W_target is given. Other components
%   are learnt in an unsupervised fashion, specifying the number of free
%   additional basis in W.
%   The matrix W_target can be estimated either using the function apl_unsupervised_nmf(),
%   storing its result W as a .mat file; or by loading an existing instrument pitch 
%   model stored in HDF5 format using the function apl_load_model().
%
% inputs: 
% - audiofile: input audio filename in WAV mono format.
% - out_folder: output filder where output audio files 
%               (*_target.wav and *_other.wav) are stored
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

%-------------------------------------------------------------
% Supervised NMF
% Load a priori basis NMF model and configuration parameters
use_piano_model = false; % set to true for assignment 3

if use_piano_model,
    % Piano model - standard NMF model (CNMF 1-column)
    nmf_model = './models/011PFNOM_norm3_1col.cnmfModel'; % NMF Model file in HDF5 format
    [W, param] = apl_load_model(nmf_model, param); % load supervised NMF model (e.g. piano model)
else
    % generic W matrix stored in apl_estimation_nmf
    load('./my_W_basis.mat'); % get previsouly estimated W basis
end

% Unsupervised basis W_other
param.n_other = 50; % number of additional "other" unsupervised basis; default = 0,



%-------------------------------------------------------------
% Informed Source Separation (optional)
% Additional constrains for NMF estimation step

% H0 constraints for activations of target and other basis (e.g. for score-informed case)
H0mask = ones(size(W,2),nFrames); % set to 0 activations for W(target) that will be inactive. (default: 1)
H0_othermask = ones(param.n_other,size(H0mask,2)); % set to 0 activations for W(other) that will be inactive (default: 1)

% generate random initialization for target and other activation gains
nBasisConstrain = size(W,2);
if param.n_other>0
    H0init = [rand(size(H0mask)) .* H0mask ; rand(param.n_other,size(H0mask,2)) .* H0_othermask];
else
    H0init = [rand(size(H0mask)) .* H0mask];
end


%-------------------------------------------------------------
% Factorization

fprintf('*** Performing factorization ***\n');

% Use W basis as fixed basis (W)
[Wout,Hout,errs,vout] = nmf_kl_con(Y',nBasisConstrain, 'niter', param.niter, 'win', param.width,'W',W,'H0',H0init,'norm_w',param.norm_w, 'n_other',param.n_other);
numBasisOut = size(Wout,2);

% target basis/activations are Wt and Ht
Wt = Wout(:,1:nBasisConstrain,:);
Ht = Hout(1:nBasisConstrain,:);
% other basis/activations initilized to zero when not used
Wo = Wout(:,nBasisConstrain+1:numBasisOut,:);
Ho = Hout(nBasisConstrain+1:numBasisOut,:);

% reconstruct estimated target and other spectrogram
% convolutive NMF
target =  rec_cnmf(Wt,Ht,param.eps); % Wt * Ht;
other = rec_cnmf(Wo,Ho,param.eps); %Wo * Ho;

% compute filter to TF input signal (Wiener filter)
maskTarget = double(target) ./ (double(target) + double(other) + eps); % avoid division by zero
maskOther = double(other) ./ (double(target) + double(other) + eps);



                                           
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

                                           
                                           
         
