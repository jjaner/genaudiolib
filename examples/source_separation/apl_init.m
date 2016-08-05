function param = apl_init()
% function param = apl_init()
% APL INITIALIZATION FUNCTION 
%
% It returns the default parameters for running the APL source separation scripts.
% 
% MTG-UPF, March 2015, Jordi Janer


% include paths
addpath('../src'); % genaudiolib
addpath('./3rdparty/nmflib');
addpath('./3rdparty/midi_lib');



%--------------------------------------------------------------------------
% CONFIGURATION PARAMETERS
% TF analysis
param.fftsize = 4096;
param.winsize = 4096;
param.hopsize = 512;%
param.windowtype = 'blackmanharris';
param.overlapWinSize = 1024;
param.signalRepresentation = 'stft';
param.zerophase = 1;

% NMF decomposition parameters
param.niter =  30; % %number of iterations , default: 30;
param.eps = 1e-20;
param.update_W = false; % set to true if target basis in W are updated. (default: false)

% MIDI file parameters
param.midi_tpqn = 120;
param.midi_bpm = 120;  

% Resynthesis parameters
filterTypes = {'Wiener', 'Sqrt-Wiener', 'Mixed'};
param.filterType = filterTypes{1}; % 3

