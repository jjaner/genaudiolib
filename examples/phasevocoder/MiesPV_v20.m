%-------- MIESPV CODE ------------
%
%  MIESPV - Monophonic audIo procESsing based on Phase Vocoder
%  This software implements a set of sound transformations based on phase
%  vocoder. The available sound transformations are adaptations of extiting implementations in order to meet the goals
%  of the sound installation authored by Edwin van der Heide in the
%  context os SONAR 2014 in Barcelona. 
%
% The MIESPV code consists of three components:
% - MiesPV_moprph.m: A Matlab script adapted to produce the sound
% transformations required for the Edwin van der Heide's installation. It loads a score
% file, and outputs a number of transformed audio files for all harmonic
% partials. This script is developed by Jordi Janer.
% - GenaudioLib library: a Matlab code library (GenAudioLib) for STFT analysis/synthesis and
% Phase Vocoding. This library is proprierty of the MTG-UPF and developed
% by Jordi Janer. 
% - YIN library: a third-party Matlab library (YIN) for pitch estimations. This library is NOT
% propierty of the MTG-UPF. http://audition.ens.fr/adc/
%
% 
% v1.3: Implements morphing between two sounds.
% v1.5: Implements reading parameters from score file.
% v1.6: Reducing HF artifacts on interpolating morphed sounds.
% v1.7: Smoothing crossfade transition for morphed outputs A/B and timbre
% correction. Pitch shift is specified in semitones.
% v1.8: Optimizations in the morphing step
% v1.9: Extra interpolation for peak selection in phase vocoder algorithm.
% v2.0: Block processing to allow long files.
%
% Collaboration of Music Technology Group/Universitat Pompeu Fabra, Sonar+D
% 2014 and Edwin van der Heide.
%
% SMC student: Antonio Sá
% MTG supervisor: Jordi Janer
% Barcelona May 2014.
% 

%  v2.0 Jordi Janer, MTG-UPF 2014
% 10/11/2014: Revision for Block-wise process to manipulate large input or output
% files without memory limitations


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER PARAMETERS

%
% input files
inputfileA = '../edwin/test_blocks/iii.wav';
inputfileB = '../edwin/test_blocks/aaa.wav';


% score file
scorefile = '../edwin/test_blocks/test_morph.score';

% output folder
outfolder = '../edwin/test_blocks/out/';
% harmonic parameters
nHarmonics = 20;%40; % 40 % number of output harmonic partials (one file per harmonic will be created)
testDuration = -1;  % limits the duration of the output sound to the specified value in seconds. Default: -1 (original size is used)
doMorphing = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INIT PARAMETERS
addpath('./yin');
addpath('../src'); % release genaudiolib

% define output filename
[p na e] = fileparts(inputfileA);
[p nb e] = fileparts(inputfileB);
[p ns e] = fileparts(scorefile);
if ~exist(outfolder)
   mkdir (outfolder);
end
if ~exist(outfolder)
mkdir (fullfile(outfolder,'temp'));
end
outfile = fullfile(outfolder,'temp',[na,'-',nb,'_',ns,'.wav']);

% Spectral analysis parameters
param.fftsize = 4096;
param.winsize = 4096;
param.hopsize = 512; %512;%
param.windowtype = 'blackmanharris';
param.zerophase = 1;
param.minPitchHz = 35;
param.maxPitchHz = 500;
param.numMaxHarmonics = 120; % for a speech pitch of 100Hz it covers 12k
useTimbrePreserve = false; % Use timbre (formant) perservation. Set to true for large pitch shift values. Default: false;

% block processing params
block_duration = 2; % 10 seconds
block_overlap = 0.25;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCORE READING
disp('reading input score...')
if exist(scorefile,'file')
   [ ~, param.fs] = wavread(inputfileA);
[synPSA, synPSB, synTimeIdxA, synTimeIdxB, morphAB] = read_score_edwin(scorefile,param);
if testDuration > -1,    
    testFr = floor(testDuration*param.fs / param.hopsize);
    synPSA = synPSA(1:testFr);
    synPSB = synPSB(1:testFr);
    synTimeIdxA = synTimeIdxA (1:testFr);
    synTimeIdxB = synTimeIdxB (1:testFr);
    morphAB = morphAB(1:testFr);
end
else
    disp('Score file not found. Exit');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get data
%-------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET DATA SOUND A

% get STFT
[Y Yphase param] = get_stft(inputfileA, param);
nFrames = size(Y,1);
% get input audio pitch curve (on OSX only works if .pitch file is available)
pitchFilename = strrep(inputfileA,'.wav','.pitch');
if ~exist(pitchFilename,'file')
    disp(['Pitch file (',pitchFilename ,') could not be found. Pitch is estimated with YIN algorithm. ']);  
 %   disp('MiesPV OSX warning: please locate the .pitch file in the same folder as the audio file. YIN is not supported for latest OSX versions.')
else
    disp(['Pitch file (',pitchFilename ,') is loaded correctly.']);  
end
pitch = get_pitch(pitchFilename, nFrames, param.hopsize, param.fs);
pitch = pitch .* (pitch < param.maxPitchHz);
% resize input data to same duration
nFrames = min(nFrames, size(pitch,1));
if  testDuration > -1 && floor(testDuration*param.fs / param.hopsize ) < nFrames,
    nFrames = floor(testDuration*param.fs / param.hopsize); % 1 second limit for generating test sounds
end

Y = Y(1:nFrames,:);
Yphase = Yphase(1:nFrames,:);
pitch = pitch(1:nFrames);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET DATA SOUND B

% get STFT
[YB YphaseB param] = get_stft(inputfileB, param);
nFramesB = size(YB,1);

% get input audio pitch curve (on OSX only works if .pitch file is available)
pitchFilenameB = strrep(inputfileB,'.wav','.pitch');
if ~exist(pitchFilenameB,'file')
    disp(['Pitch file (',pitchFilenameB ,') could not be found. ']);  
    disp('MiesPV OSX warning: please locate the .pitch file in the same folder as the audio file. YIN is not supported for latest OSX versions.')
else
    disp(['Pitch file (',pitchFilenameB ,') is loaded correctly.']);  
end
pitchB = get_pitch(pitchFilenameB, nFramesB, param.hopsize, param.fs);
pitchB = pitchB .* (pitchB < param.maxPitchHz);

% resize input data to same duration
nFramesB = min(nFramesB, size(pitchB,1));
if  testDuration > -1 && floor(testDuration*param.fs / param.hopsize ) < nFramesB,
    nFramesB = floor(testDuration*param.fs / param.hopsize); % 1 second limit for generating test sounds
end

YB = YB(1:nFramesB,:);
YphaseB = YphaseB(1:nFramesB,:);
pitchB = pitchB(1:nFramesB);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute morphing data
%-------------------------------------------------------------
% compute target transformation parameters
nFramesC = min(size(synTimeIdxA,1),size(synTimeIdxB,1));
morphFactor = 1 - morphAB(1:nFramesC); % inverse so that morphAB=0 means all to A
itype =  'linear'; % 'spline'; 
pitchAinterp = interp1([1:size(pitch,1) ],pitch,synTimeIdxA(1:nFramesC),itype);
pitchBinterp = interp1([1:size(pitchB,1) ],pitchB,synTimeIdxB(1:nFramesC),itype);
% synpitchA  = synPSA(1:nFramesC)' .* (pitch(round(synTimeIdxA(1:nFramesC))));
% synpitchB  = synPSB(1:nFramesC)' .* (pitchB(round(synTimeIdxB(1:nFramesC))));
synpitchA  = synPSA(1:nFramesC)' .* pitchAinterp;
synpitchB  = synPSB(1:nFramesC)' .* pitchBinterp; 

% target pitch 
synpitchA_cont = get_continuous_pitch(synpitchA);
synpitchB_cont = get_continuous_pitch(synpitchB);

synPitchC  = morphFactor .* synpitchA_cont + (1- morphFactor).* synpitchB_cont;
% synPSAC = (synPitchC./pitch(round(synTimeIdxA(1:nFramesC)))) .* (synpitchA>0) + (synpitchA<=0); % set PS to 1 if pitch == 0
% synPSBC = (synPitchC./pitchB(round(synTimeIdxB(1:nFramesC)))) .* (synpitchB>0) + (synpitchB<=0);% set PS to 1 if pitch == 0
synPSAC = ones(size(synPitchC));
validIDx = (pitchAinterp>0);
synPSAC(validIDx) = (synPitchC(validIDx)./pitchAinterp(validIDx)); % set PS to 1 if pitch == 0
synPSBC = ones(size(synPitchC));
validIDx = (pitchBinterp>0);
synPSBC(validIDx) = (synPitchC(validIDx)./pitchBinterp(validIDx));% set PS to 1 if pitch == 0

if ~doMorphing
    synPSAC = synPSA;
    synPSBC = synPSB;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iterative block processing - blocks stored
%-------------------------------------------
disp('Start block processing...');


block_duration_frames = floor(block_duration *param.fs / param.hopsize);
block_overlap_frames =  floor(block_overlap*param.fs / param.hopsize);
% round block duration in seconds to match multiple of frame duration 
block_duration2 = block_duration_frames * param.hopsize / param.fs ;

nBlocks = ceil((nFramesC ) / (block_duration_frames));
nFramesCTotal = nFramesC;


if 1 % skip process for test concatenation
    
% start block processing
for block=0:nBlocks-1,
    
    % modify variables per block
    bFr =  1 + (((block) * block_duration_frames) );
    eFr =  1 + (( ((block+1) * block_duration_frames + block_overlap_frames )));
    eFr = min(nFramesCTotal, eFr);    
    
    disp(['Processing Block: ', num2str(block)]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % each block processing
    %-------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PROCESSING SOUND A
    
    % transform object (pitch shifting and time scaling)
    [outmagA, outphaseA, outpitchA] = compute_phasevoc(Y, Yphase, param, pitch, synPSAC (bFr:eFr,:), synTimeIdxA (bFr:eFr,:), useTimbrePreserve);
    
    
    % % compute residual for A
    harmMaskResidual = createHarmonicPartialsMask(outpitchA, [0, 1, nHarmonics], param);
    [outResMagA outResPhaseA] = computeResidualFromHarmMask(Y(min(round(synTimeIdxA (bFr:eFr,:)),nFrames),:), ...
        outmagA, outphaseA, pitch(min(round(synTimeIdxA (bFr:eFr,:)),nFrames)), harmMaskResidual, nHarmonics, param);
    
    disp('Sound A processed.');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PROCESSING SOUND B
    
    % transform object (pitch shifting and time scaling)
    [outmagB, outphaseB, outpitchB] = compute_phasevoc(YB, YphaseB, param, pitchB, synPSBC (bFr:eFr,:), synTimeIdxB (bFr:eFr,:), useTimbrePreserve);
    
    % % compute residual for B
    harmMaskResidual = createHarmonicPartialsMask(outpitchB, [0, 1, nHarmonics], param);
    [outResMagB outResPhaseB] = computeResidualFromHarmMask(YB(min(round(synTimeIdxB (bFr:eFr,:)),nFramesB),:), ...
        outmagB, outphaseB, pitchB(min(round(synTimeIdxB (bFr:eFr,:)),nFramesB)), harmMaskResidual, nHarmonics, param);
    
    disp('Sound B processed.');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PROCESSING MORPH A-B
    
    nFramesC = size(bFr:eFr,2); % new v20. Test
    
    durFramesC = nFramesC;
    outmagA = outmagA(1:nFramesC,:);
    outphaseA = outphaseA(1:nFramesC,:);
    outmagB = outmagB(1:nFramesC,:);
    outphaseB = outphaseB(1:nFramesC,:);
    energyA = sum(outmagA.^2,2);
    energyB = sum(outmagB.^2,2);
    % save outputs before morphing
    outfileC = strrep(outfile,'.wav',['_fromA','_block_' num2str(block),'.wav']);
    write_stft(outmagA, outphaseA, outfileC, param);
    outfileC = strrep(outfile,'.wav',['_fromB','_block_' num2str(block),'.wav']);
    write_stft(outmagB, outphaseB, outfileC, param);
    
    % target timbre
    [outmagEnv_A, ~, ~, ~] = harmonicEnv(outmagA, outpitchA, param.fs, param.fs/2, param.numMaxHarmonics);
    [outmagEnv_B, ~, ~, ~] = harmonicEnv(outmagB, outpitchB, param.fs, param.fs/2, param.numMaxHarmonics);
    timberC = ((morphFactor(bFr:eFr)./sqrt(energyA)) * ones(1,size(outmagEnv_A,2))) .* (outmagEnv_A(1:nFramesC,:)) +  (((1- morphFactor (bFr:eFr) )./sqrt(energyB)) * ones(1,size(outmagEnv_A,2))) .* outmagEnv_B(1:nFramesC,:);
    
    
    % Transform object (pitch shifting and time scaling)
    % process output A
    envFilter = computeTimbreFilter(timberC, outmagA, outpitchA,  param);
    outmagC = outmagA .* envFilter;
    outphaseC = outphaseA;
    % process output B
    envFilter2 = computeTimbreFilter(timberC, outmagB, outpitchB,  param);
    outmagC2 = outmagB  .* (envFilter2);
    outphaseC2 = outphaseB;
    clearvars 'outmagA outmagB outphaseA outphaseB';
    
    % save morphed signals
    outfileC = strrep(outfile,'.wav',['_morphfromA','_block_' num2str(block),'.wav']);
    write_stft(outmagC, outphaseC, outfileC, param);
    outfileC = strrep(outfile,'.wav',['_morphfromB','_block_' num2str(block),'.wav']);
    write_stft(outmagC2, outphaseC2, outfileC, param);
    
    % output mixture of sounds / cross-fade for transition
    outmag = zeros(size(outmagC));
    outphase = zeros(size(outmagC));
    crossFadeMorph = sigmf(morphFactor (bFr:eFr) ,[20 0.5]); % non-linear (signmoid) sound source selection (A, B) for better sound consistency at the extreme morphing values.
    for i = 1:size(crossFadeMorph,1),
        % complex crossfade
        [outmag(i,:), outphase(i,:)] = c_spectral_mix(outmagC(i,:), outmagC2(i,:), outphaseC(i,:), outphaseC2(i,:), crossFadeMorph(i), (1-crossFadeMorph(i)) );
        [outResMag(i,:), outResPhase(i,:)] = c_spectral_mix(outResMagA(i,:), outResMagB(i,:), outResPhaseA(i,:), outResPhaseB(i,:), crossFadeMorph(i), (1-crossFadeMorph(i)) );
        
    end
    outfileC = strrep(outfile,'.wav',['_morph_mix','_block_' num2str(block),'.wav']);
    write_stft(outmag, outphase, outfileC, param);
    clearvars 'outmagC outmagC2 outphaseC outphaseC2';
    
    disp('Morphing outputs created.')
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT HARMONIC PARTIALS DECOMPOSITON
    disp('Computing harmonic partial outputs...')
    % partial filtering process
    harmonicMask = zeros(size(outmag));
    outpitch = synPitchC (bFr:eFr) .* (outpitchA>0 | outpitchB>0); % Output pitch only if any of both
    
    % Geneate filtered outputs one file per harmonic
    for hh=1:nHarmonics % 20 harmonic partials
        
        % create time-varying mask to filter only certain patrials
        % harmonics_list variable contains list of triples for by-passing partials:
        % row: {timeInFrames, firstPartial, lastPartial}
        harmonics_list = [ ...
            0*nFramesC , hh , hh; % filter single partial during the whole file duration
            ];
        
        partialsMask = createHarmonicPartialsMask(outpitch, harmonics_list, param);
        outmag2 = outmag .* partialsMask;
        harmonicMask = min(1,(harmonicMask + partialsMask));
        
        % write output audio
        outfile2 = strrep(outfile,'.wav',['_',num2str(hh),'_block_' num2str(block),'.wav']);
        write_stft(outmag2, outphase, outfile2, param)
        
        
    end % loop generate output per partial
    
    % computed non-harmonic components: "rest" and "unvoiced"
    restMask = (1 - harmonicMask);
    unvoicedMask = (1 - harmonicMask) .*((outpitch<param.minPitchHz) * ones(1,size(harmonicMask,2)));
    outmag3 = outmag.*restMask;
    outmag4 = outmag .*unvoicedMask;
    outmag5 = outmag .*harmonicMask;
    
    % save inverse mask (rest)
    outfile3 = strrep(outfile,'.wav',['_invharmonic','_block_' num2str(block),'.wav']);
    write_stft(outmag3, outphase, outfile3, param)
    % save unvoiced frames
    outfile4 = strrep(outfile,'.wav',['_unvoiced','_block_' num2str(block),'.wav']);
    write_stft(outmag4, outphase, outfile4, param)
    % save harmonic frames
    outfile5 = strrep(outfile,'.wav',['_harmonic','_block_' num2str(block),'.wav']);
    write_stft(outmag5, outphase, outfile5, param)
    % save total output
    outfile6 = strrep(outfile,'.wav',['_total','_block_' num2str(block),'.wav']);
    write_stft(outmag, outphase, outfile6, param)
    % estimate residual
    outfile7 = strrep(outfile,'.wav',['_residual','_block_' num2str(block),'.wav']);
    write_stft(outResMag, outResPhase, outfile7, param);
    
    
end % loop for block processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%clear all;
end % skip process test concat 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concatenate block files
disp('Concatenate blocks and write output files...')

% save mix signal
outfile = fullfile(outfolder,[na,'-',nb,'_',ns,'.wav']);
outfolder = fullfile(outfolder,'temp');

outfileC = strrep(outfile,'.wav',['_morph_mix','.wav']);
[p sourceName e] = fileparts(outfileC);
concatenateBlocks(outfolder, outfileC, sourceName, param.hopsize, block_duration2, block_overlap)

    
% save morphed signals
outfileC = strrep(outfile,'.wav',['_morphfromA','.wav']);
[p sourceName e] = fileparts(outfileC);
concatenateBlocks(outfolder, outfileC, sourceName, param.hopsize, block_duration2, block_overlap)


outfileC = strrep(outfile,'.wav',['_morphfromB','.wav']);
[p sourceName e] = fileparts(outfileC);
concatenateBlocks(outfolder, outfileC, sourceName, param.hopsize, block_duration2, block_overlap)


% harmonic partials
for hh=1:nHarmonics % 20
    outfileC = strrep(outfile,'.wav',['_',num2str(hh),'.wav']);
    [p sourceName e] = fileparts(outfileC);
    concatenateBlocks(outfolder, outfileC, sourceName, param.hopsize, block_duration2, block_overlap)
end

% save inverse mask (rest)
outfileC = strrep(outfile,'.wav',['_invharmonic','.wav']);
[p sourceName e] = fileparts(outfileC);
concatenateBlocks(outfolder, outfileC, sourceName, param.hopsize, block_duration2, block_overlap)


% save unvoiced frames
outfileC = strrep(outfile,'.wav',['_unvoiced','.wav']);
[p sourceName e] = fileparts(outfileC);
concatenateBlocks(outfolder, outfileC, sourceName, param.hopsize, block_duration2, block_overlap)

% save harmonic frames
outfileC = strrep(outfile,'.wav',['_harmonic','.wav']);
[p sourceName e] = fileparts(outfileC);
concatenateBlocks(outfolder, outfileC, sourceName, param.hopsize, block_duration2, block_overlap)

% save total output
outfileC = strrep(outfile,'.wav',['_total','.wav']);
[p sourceName e] = fileparts(outfileC);
concatenateBlocks(outfolder, outfileC, sourceName, param.hopsize, block_duration2, block_overlap)

% estimate residual
outfileC = strrep(outfile,'.wav',['_residual','.wav']);
[p sourceName e] = fileparts(outfileC);
concatenateBlocks(outfolder, outfileC, sourceName, param.hopsize, block_duration2, block_overlap)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
